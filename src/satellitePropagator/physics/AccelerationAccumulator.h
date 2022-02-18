//
// Created by Oliver on 12.05.21.
//

#pragma once
#include <cmath>

#include <array>
#include <numeric>

#include "Constants.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/io/FileOutput.h"
#include "satellitePropagator/physics/AccelerationComponents/include.h"

/**
 * @namespace Acceleration
 *
 * @brief Contains classes and functions to calculate accelerations for
 * Debris::Debris objects
 *
 * The implementations are following the equations in <a href="Math.pdf"> math reference</a>
 *
 */
namespace Acceleration {
/**
 * @class AccelerationAccumulator
 *
 * @brief Accumulates all accelerations of the
 * Acceleration::AccelerationComponent specified by the #config vector
 */
template <class Container>
class AccelerationAccumulator {
public:
    /**
     * @brief Default constructor
     *
     * creates ann AccelerationAccumulator object with all values initialized to zero
     */
    AccelerationAccumulator();
    /**
     * @brief Constructor creates a new Acceleration::AcceleratorAccumulator
     * object
     *
     * @param config_arg The 8D bool vector encoding the
     * Acceleration::AccelerationComponent to apply in the simulation. The Order
     * of flags is [::KEP, ::J2, ::C22, ::S22, ::SOL, ::LUN, ::SRP, ::DRAG]
     * @param container Reference to the Container object holding
     * the Debris::Debris objects to apply acceleration to
     * @param t_arg Current time
     * @param file_output_arg Reference to the FileOutput object used for output of the acceleration values during computation
     */
    AccelerationAccumulator(const std::array<bool, 8>& config_arg,
        Container& container,
        double t_arg,
        FileOutput<Container>* file_output_arg)
        : config(config_arg)
        , container(&container)
        , t(t_arg)
        , file_output(file_output_arg) {};

    /**
     * @brief Default destructor
     *
     * Destroys the Acceleration::AccelerationAccumulator object
     */
    ~AccelerationAccumulator();

    /**
     * @brief Applies all Acceleration::AccelerationComponent specified by
     * #config
     *
     * Loops over all Debris::Debris objects in #container.
     * Checks if the corresponding flag in the #config vector is set to true and
     * calls the apply() functions. Sets the #Debris::Debris::acc_t1 value of
     * every Debris::Debris object.
     *
     * The possible apply() functions are:
     * - Acceleration::KepComponent::apply() if #config[Acceleration::KEP]
     * - Acceleration::J2Component::apply() if #config[Acceleration::J2]
     * - Acceleration::C22Component::apply() if #config[Acceleration::C22]
     * - Acceleration::S22Component::apply() if #config[Acceleration::S22]
     * - Acceleration::SolComponent::apply() if #config[Acceleration::SOL]
     * - Acceleration::LunComponent::apply() if #config[Acceleration::LUN]
     * - Acceleration::SRPComponent::apply() if #config[Acceleration::SRP]
     * - Acceleration::DragComponent::apply() if #config[Acceleration::DRAG]
     *
     * @tparam write If true triggers this class' file output writer.
     */
    template <bool write>
    void applyComponents() const;

    /**
     * @brief does the same as #applyComponents() plus output
     *
     * Calculates the needed acceleration components and writes the value for each one to a csv file
     */
    void applyAmdWriteComponents() const;

private:
    /**
     * @brief 8D bool vector encoding the Acceleration::AccelerationComponent to
     * apply in the simulation.
     *
     * 8D bool vector encoding the Acceleration::AccelerationComponent to apply
     * in the simulation. The Order of flags is
     * - #Acceleration::KEP
     * - #Acceleration::J2
     * - #Acceleration::C22
     * - #Acceleration::S22
     * - #Acceleration::SOL
     * - #Acceleration::LUN
     * - #Acceleration::SRP
     * - #Acceleration::DRAG
     */
    std::array<bool, 8> config {};
    Container*
        container
        = nullptr; /**< Reference to the Container object holding the
             Debris::Debris objects to apply acceleration to*/
    double t = 0; /**< current time*/
    FileOutput<Container>* file_output = nullptr; /**< used to write detailed output data during calculations */
public:
    /**
     * @brief Getter function for #config
     *
     * @return Value of #config
     */
    [[nodiscard]] const std::array<bool, 8>& getConfig() const;
    std::array<bool, 8>& getConfig();

    /**
     * @brief Setter function for #config
     *
     * @param config New value of #config
     */
    void setConfig(const std::array<bool, 8>& config);

    /**
     * @brief Getter function for #container
     *
     * @return Value of #container
     */
    [[nodiscard]] const Container& getContainer() const;
    Container& getContainer();

    /**
     * @brief Setter function for #container
     *
     * @param container New value of #container
     */
    void setContainer(Container& container);

    /**
     * @brief Getter function for #t
     *
     * @return Value of #t
     */
    [[nodiscard]] double getT() const;

    /**
     * @brief Setter function for #t
     *
     * @param t New value of #t
     */
    void setT(double t);

    /**
     * @brief Getter function for #file_output
     *
     * @return Value of #file_output
     */
    [[nodiscard]] const FileOutput<Container>& getFileOutput() const;
    FileOutput<Container>& getFileOutput();

    /**
     * @brief Setter function for #file_output
     *
     * @param fileOutput New value of #file_output
     */
    void setFileOutput(FileOutput<Container>& fileOutput);
};
} // namespace Acceleration

//
// Created by Oliver on 12.05.21.
// based on https://www.overleaf.com/project/608918003323352d4d3ceb55
// equations. for reference
namespace Acceleration {

template <class Container>
AccelerationAccumulator<Container>::AccelerationAccumulator() = default;

template <class Container>
AccelerationAccumulator<Container>::~AccelerationAccumulator() = default;

template <class Container>
template <bool write>
void AccelerationAccumulator<Container>::applyComponents() const
{
#ifdef AUTOPAS_OPENMP
    // autopas works with parallel iterators which it controls itself. No pragma omp for needed!
#pragma omp parallel
#endif
    {
        // will be modified by the apply functions
        std::array<double, 3> new_acc_total { 0, 0, 0 };
        std::array<double, 3> new_acc_component { 0, 0, 0 };
        double d_srp { 0. };
        // are constant for this time step
        double c_term;
        double s_term;
        std::array<double, 6> sun_params {};
        // setup only needed for SolComponent and SRPComponent
        if (config[SOL] || config[SRP]) {
            sun_params = SolComponent::setUp(t);
        }
        std::array<double, 6> moon_params {};
        // setup only needed for LunComponent
        if (config[LUN]) {
            moon_params = LunComponent::setUp(t);
        }
        // setup only needed for C22Component and S22Component
        if (config[C22] || config[S22]) {
            // Eq 15
            c_term = std::cos((Physics::THETA_G + Physics::NU_EARTH * t) * Physics::RAD_FACTOR);
            s_term = std::sin((Physics::THETA_G + Physics::NU_EARTH * t) * Physics::RAD_FACTOR);
        }

        for (auto& d : *container) {
            if (write) {
                file_output->writeAcc_start(t);
            }
            new_acc_total[0] = 0;
            new_acc_total[1] = 0;
            new_acc_total[2] = 0;
            d_srp = 0;
            // Eq 1
            if (config[KEP]) {
                new_acc_component = KepComponent::apply(d);
                if (write) {
                    file_output->writeAcc_value(new_acc_component);
                }
                new_acc_total[0] += new_acc_component[0];
                new_acc_total[1] += new_acc_component[1];
                new_acc_total[2] += new_acc_component[2];
            }
            // Managed particles compensate for all perturbations and stay on the keplerian orbit. We simulate this by not applying any pertubations on them.
            if(d.getActivityState() != Container::Particle_t::ActivityState::evasivePreserving) {
                if (config[J2]) {
                    new_acc_component = J2Component::apply(d);
                    if (write) {
                        file_output->writeAcc_value(new_acc_component);
                    }
                    new_acc_total[0] += new_acc_component[0];
                    new_acc_total[1] += new_acc_component[1];
                    new_acc_total[2] += new_acc_component[2];
                }
                // if we want to calculate both C22 and S22 we can share many of the
                // calculation steps
                if (config[C22] && config[S22]) {
                    new_acc_component = C22S22Component::apply(d, c_term, s_term);
                    if (write) {
                        file_output->writeAcc_value(new_acc_component);
                        file_output->writeAcc_value(new_acc_component);
                    }
                    new_acc_total[0] += new_acc_component[0];
                    new_acc_total[1] += new_acc_component[1];
                    new_acc_total[2] += new_acc_component[2];
                } // if only one of the two should be calculated we only calculate the
                // needed one
                else {
                    if (config[C22]) {
                        new_acc_component = C22Component::apply(d, c_term, s_term);
                        if (write) {
                            file_output->writeAcc_value(new_acc_component);
                        }
                        new_acc_total[0] += new_acc_component[0];
                        new_acc_total[1] += new_acc_component[1];
                        new_acc_total[2] += new_acc_component[2];
                    }
                    if (config[S22]) {
                        new_acc_component = S22Component::apply(d, c_term, s_term);
                        if (write) {
                            file_output->writeAcc_value(new_acc_component);
                        }
                        new_acc_total[0] += new_acc_component[0];
                        new_acc_total[1] += new_acc_component[1];
                        new_acc_total[2] += new_acc_component[2];
                    }
                }
                if (config[SOL]) {
                    new_acc_component = SolComponent::apply(d, d_srp, sun_params);
                    if (write) {
                        file_output->writeAcc_value(new_acc_component);
                    }
                    new_acc_total[0] += new_acc_component[0];
                    new_acc_total[1] += new_acc_component[1];
                    new_acc_total[2] += new_acc_component[2];
                }
                if (config[LUN]) {
                    new_acc_component = LunComponent::apply(d, moon_params);
                    if (write) {
                        file_output->writeAcc_value(new_acc_component);
                    }
                    new_acc_total[0] += new_acc_component[0];
                    new_acc_total[1] += new_acc_component[1];
                    new_acc_total[2] += new_acc_component[2];
                }
                if (config[SRP]) {
                    new_acc_component = SRPComponent::apply(d, d_srp, sun_params);
                    if (write) {
                        file_output->writeAcc_value(new_acc_component);
                    }
                    new_acc_total[0] += new_acc_component[0];
                    new_acc_total[1] += new_acc_component[1];
                    new_acc_total[2] += new_acc_component[2];
                }
                if (config[DRAG]) {
                    new_acc_component = DragComponent::apply(d);
                    if (write) {
                        file_output->writeAcc_value(new_acc_component);
                    }
                    new_acc_total[0] += new_acc_component[0];
                    new_acc_total[1] += new_acc_component[1];
                    new_acc_total[2] += new_acc_component[2];
                }
            }
            d.setAccT1(new_acc_total);
            if (write) {
                file_output->writeAcc_end(new_acc_total);
            }
        }
    }
}

template <class Container>
const std::array<bool, 8>& AccelerationAccumulator<Container>::getConfig() const
{
    return config;
}

template <class Container>
std::array<bool, 8>& AccelerationAccumulator<Container>::getConfig()
{
    return config;
}
template <class Container>
void AccelerationAccumulator<Container>::setConfig(const std::array<bool, 8>& config)
{
    AccelerationAccumulator<Container>::config = config;
}

template <class Container>
const Container& AccelerationAccumulator<Container>::getContainer() const
{
    return *container;
}

template <class Container>
Container& AccelerationAccumulator<Container>::getContainer()
{
    return *container;
}
template <class Container>
void AccelerationAccumulator<Container>::setContainer(Container& container)
{
    AccelerationAccumulator<Container>::container = &container;
}

template <class Container>
double AccelerationAccumulator<Container>::getT() const
{
    return t;
}

template <class Container>
void AccelerationAccumulator<Container>::setT(double t)
{
    AccelerationAccumulator<Container>::t = t;
}

template <class Container>
const FileOutput<Container>& AccelerationAccumulator<Container>::getFileOutput() const
{
    return *file_output;
}

template <class Container>
FileOutput<Container>& AccelerationAccumulator<Container>::getFileOutput()
{
    return *file_output;
}

template <class Container>
void AccelerationAccumulator<Container>::setFileOutput(FileOutput<Container>& fileOutput)
{
    file_output = &fileOutput;
}

} // namespace Acceleration

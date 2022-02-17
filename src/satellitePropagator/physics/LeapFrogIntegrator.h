//
// Created by Oliver on 13.05.21.
//

#pragma once
#include "AccelerationAccumulator.h"
#include "satellitePropagator/debris/AccelerationUpdate.h"

/**
 * @class Integrator
 *
 * @brief Calculates status of Debris::Debris objects for next time step
 */
template <class Container>
class LeapFrogIntegrator {
public:
    /**
     * @brief Default constructor
     *
     * Creates a new Integrator object with all values zero initialized
     */
    LeapFrogIntegrator();

    /**
     * @brief Creates a Integrator object and sets the private #container and
     * #delta_t member variables
     *
     * @param container Reference to the Container object holding
     * the Debris::Debris objects to integrate for
     * @param accumulator_arg Reference to the
     * Acceleration::AccelerationAccumulator object to calculate acceleration for
     * the current time step
     * @param delta_t_arg Time step to Integrate over
     */
    LeapFrogIntegrator(Container& container,
        Acceleration::AccelerationAccumulator<Container>& accumulator_arg,
        double delta_t_arg)
        : container(&container)
        , accumulator(&accumulator_arg)
        , delta_t(delta_t_arg) {};

    /**
     * @brief Default destructor
     *
     * Destroys the Integrator object
     */
    virtual ~LeapFrogIntegrator();

    /**
     * @brief Calculates a complete time step
     *
     * Integrator::calculateAcceleration() is called first.
     * Because the Integrator::calculatePosition() function needs the
     * Debris::Debris::velocity vector of the last time step of the
     * Debris::Debris objects of the Container #container it is called
     * before Integrator::calculateVelocity()
     *
     * @param write_time_step if true all calculated acceleration components are writen to a file
     */
    void integrate(bool write_time_step = false) const;

    /**
     * @brief Calculates the new position
     *
     * Calculates the Debris::Debris::position vector for all Debris::Debris
     * objects of the Container #container Uses leapfrog integration
     * with Debris::Debris::velocity and Debris::Debris::acc_t0
     *
     */
    void calculatePosition() const;

    /**
     * @brief Calculates the new velocities
     *
     * Calculates the Debris::Debris::velocity vector for all Debris::Debris
     * objects of the Container #container Uses leapfrog integration
     * with Debris::Debris::acc_t0 and Debris::Debris::acc_t1
     *
     */
    void calculateVelocity() const;

    /**
     * @brief Calculates te acceleration for the current time step
     *
     * Calls the Acceleration::AccelerationAccumulator::applyComponents()
     * function of the Acceleration::AccelerationAccumulator object #accumulator
     * After this function is called the Debris::Debris objects of the
     * Container #container are ready for integrating
     *
     * @param write_time_step if true all calculated acceleration components are writen to a file
     */
    void calculateAcceleration(bool write_time_step) const;

private:
    Container*
        container
        = nullptr; /**< Reference to the Container object holding the
             Debris::Debris objects to integrate for*/
    Acceleration::AccelerationAccumulator<Container>*
        accumulator
        = nullptr; /**< Reference to the Acceleration::AccelerationAccumulator
                   object to calculate acceleration for the current time
                   step*/
    double delta_t = 0; /**< Time step to Integrate over */
public:
    /**
     * @brief Getter function for #delta_t
     *
     * @return Value of #delta_t
     */
    [[nodiscard]] double getDeltaT() const;

    /**
     * @brief Setter function for #delta_t
     *
     * @param deltaT New value of #delta_t
     */
    void setDeltaT(double deltaT);

    /**
     * @brief Getter function for #accumulator
     *
     * @return Value of #accumulator
     */
    [[nodiscard]] const Acceleration::AccelerationAccumulator<Container>& getAccumulator() const;
    Acceleration::AccelerationAccumulator<Container>& getAccumulator();

    /**
     * @brief Setter function for #accumulator
     *
     * @param accumulator New value of #accumulator
     */
    void setAccumulator(Acceleration::AccelerationAccumulator<Container>& accumulator);

    /**
     * @brief Getter function for #container
     *
     * @return Value of #container
     */
    [[nodiscard]] const Container& getContainer() const;
    Container& getContainer();

    /**
     * @brief Setter function for #delta_t
     *
     * @param container New value of #delta_t
     */
    void setContainer(Container& container);
};

template <class Container>
LeapFrogIntegrator<Container>::LeapFrogIntegrator() = default;

template <class Container>
LeapFrogIntegrator<Container>::~LeapFrogIntegrator() = default;

template <class Container>
void LeapFrogIntegrator<Container>::integrate(bool write_time_step) const
{
    AccelerationUpdate::accelerationUpdate(container);
    calculatePosition();
    accumulator->setT(accumulator->getT() + delta_t);
    calculateAcceleration(write_time_step);
    calculateVelocity();
    // update time
}

template <class Container>
void LeapFrogIntegrator<Container>::calculatePosition() const
{
#ifdef AUTOPAS_OPENMP
    // autopas works with parallel iterators which it controls itself. No pragma omp for needed!
#pragma omp parallel
#endif
    {
        double factor = delta_t * delta_t * 0.5;
        std::array<double, 3> new_pos {};
        for (auto& d : *container) {
            new_pos = d.getPosition();
            new_pos[0] = new_pos[0] + delta_t * d.getVelocity()[0] + factor * d.getAccT0()[0];
            new_pos[1] = new_pos[1] + delta_t * d.getVelocity()[1] + factor * d.getAccT0()[1];
            new_pos[2] = new_pos[2] + delta_t * d.getVelocity()[2] + factor * d.getAccT0()[2];
            d.setPosition(new_pos);
        }
    }
}

template <class Container>
void LeapFrogIntegrator<Container>::calculateVelocity() const
{
#ifdef AUTOPAS_OPENMP
    // autopas works with parallel iterators which it controls itself. No pragma omp for needed!
#pragma omp parallel
#endif
    {
        const double factor = delta_t * 0.5;
        std::array<double, 3> new_velocity {};
        for (auto& d : *container) {
            new_velocity = d.getVelocity();
            new_velocity[0] = new_velocity[0] + factor * (d.getAccT0()[0] + d.getAccT1()[0]);
            new_velocity[1] = new_velocity[1] + factor * (d.getAccT0()[1] + d.getAccT1()[1]);
            new_velocity[2] = new_velocity[2] + factor * (d.getAccT0()[2] + d.getAccT1()[2]);
            d.setVelocity(new_velocity);
        }
    }
}

template <class Container>
void LeapFrogIntegrator<Container>::calculateAcceleration(bool write_time_step) const
{
    if (write_time_step) {
        accumulator->template applyComponents<true>();
    } else {
        accumulator->template applyComponents<false>();
    }
}

template <class Container>
double LeapFrogIntegrator<Container>::getDeltaT() const
{
    return delta_t;
}

template <class Container>
void LeapFrogIntegrator<Container>::setDeltaT(double deltaT)
{
    delta_t = deltaT;
}

template <class Container>
const Container& LeapFrogIntegrator<Container>::getContainer() const
{
    return *container;
}

template <class Container>
Container& LeapFrogIntegrator<Container>::getContainer()
{
    return *container;
}

template <class Container>
void LeapFrogIntegrator<Container>::setContainer(Container& container)
{
    LeapFrogIntegrator<Container>::container = &container;
}

template <class Container>
const Acceleration::AccelerationAccumulator<Container>& LeapFrogIntegrator<Container>::getAccumulator() const
{
    return *accumulator;
}

template <class Container>
Acceleration::AccelerationAccumulator<Container>& LeapFrogIntegrator<Container>::getAccumulator()
{
    return *accumulator;
}

template <class Container>
void LeapFrogIntegrator<Container>::setAccumulator(Acceleration::AccelerationAccumulator<Container>& accumulator)
{
    LeapFrogIntegrator<Container>::accumulator = &accumulator;
}

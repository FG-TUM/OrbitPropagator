//
// Created by FG-TUM on 17.02.22.
//

#pragma once
#include <cmath>

#include "AccelerationAccumulator.h"
#include "satellitePropagator/debris/AccelerationUpdate.h"
/**
 * @class Integrator
 *
 * @brief Calculates status of Debris::Debris objects for next time step
 */
template <class Container>
class YoshidaIntegrator {
public:
    /**
     * @brief Default constructor
     *
     * Creates a new Integrator object with all values zero initialized
     */
    YoshidaIntegrator();

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
    YoshidaIntegrator(Container& container,
        Acceleration::AccelerationAccumulator<Container>& accumulator_arg,
        double delta_t_arg)
        : container(&container)
        , accumulator(&accumulator_arg)
        , delta_t(delta_t_arg)
        , delta_tDiv4(delta_t_arg / 4.){};

    /**
     * @brief Default destructor
     *
     * Destroys the Integrator object
     */
    virtual ~YoshidaIntegrator();

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
    void calculatePosition(double c) const;

    /**
     * @brief Calculates the new velocities
     *
     * Calculates the Debris::Debris::velocity vector for all Debris::Debris
     * objects of the Container #container Uses leapfrog integration
     * with Debris::Debris::acc_t0 and Debris::Debris::acc_t1
     *
     */
    void calculateVelocity(double d) const;

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
    double delta_t; /**< Time step to Integrate over */
    double delta_tDiv4; /**< Time step divided by 4 */

    // constants
    const std::array<double, 2> w {
        -std::cbrt(2.) / (2. - std::cbrt(2.)),
        1 / (2. - std::cbrt(2.)),
    };
    const std::array<double, 4> c {
        w[1] / 2.,
        (w[0] + w[1]) / 2.,
        (w[0] + w[1]) / 2.,
        w[1] / 2.,
    };
    const std::array<double, 3> d {
        w[1],
        w[0],
        w[1],
    };

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
YoshidaIntegrator<Container>::YoshidaIntegrator() = default;

template <class Container>
YoshidaIntegrator<Container>::~YoshidaIntegrator() = default;

template <class Container>
void YoshidaIntegrator<Container>::integrate(bool write_time_step) const
{
    AccelerationUpdate::accelerationUpdate(container);
    // do three full and one half sub step
    for (size_t i = 0; i < 3; ++i) {
        calculatePosition(c[i]);
        accumulator->setT(accumulator->getT() + delta_tDiv4);
        calculateAcceleration(write_time_step);
        calculateVelocity(d[i]);
    }
    calculatePosition(c[3]);
    accumulator->setT(accumulator->getT() + delta_tDiv4);

    // Le Math:
    // x_i^1 = x_i + c_1 v_i delta t
    // v_i^1 = v_i + d_1 a(x_i^1) delta t
    // x_i^2 = x_i^1 + c_2 v_i^1 delta t
    // v_i^2 = v_i^1 + d_2 a(x_i^2) delta t
    // x_i^3 = x_i^2 + c_3 v_i^2 delta t
    // v_i^3 = v_i^2 + d_3 a(x_i^3) delta t
    // x_i^4 = x_i^3 + c_4 v_i^3 delta t        = x_(i+1)
    // v_i^4 = v_i^3                            = v_(i+1)
}

template <class Container>
void YoshidaIntegrator<Container>::calculatePosition(double c) const
{
#ifdef AUTOPAS_OPENMP
    // autopas works with parallel iterators which it controls itself. No pragma omp for needed!
#pragma omp parallel
#endif
    {
        for (auto& particle : *container) {
            auto v = particle.getVelocity();
            auto pos = particle.getPosition();
            pos[0] += c * delta_t * v[0];
            pos[1] += c * delta_t * v[1];
            pos[2] += c * delta_t * v[2];
            particle.setPosition(pos);
        }
    }
}

template <class Container>
void YoshidaIntegrator<Container>::calculateVelocity(double d) const
{
#ifdef AUTOPAS_OPENMP
    // autopas works with parallel iterators which it controls itself. No pragma omp for needed!
#pragma omp parallel
#endif
    {
        for (auto& particle : *container) {
            auto v = particle.getVelocity();
            const auto& a = particle.getAccT0();
            v[0] += d * a[0] * delta_t;
            v[1] += d * a[1] * delta_t;
            v[2] += d * a[2] * delta_t;
            particle.setVelocity(v);
        }
    }
}

template <class Container>
void YoshidaIntegrator<Container>::calculateAcceleration(bool write_time_step) const
{
    if (write_time_step) {
        accumulator->template applyComponents<true>();
    } else {
        accumulator->template applyComponents<false>();
    }
}

template <class Container>
double YoshidaIntegrator<Container>::getDeltaT() const
{
    return delta_t;
}

template <class Container>
void YoshidaIntegrator<Container>::setDeltaT(double deltaT)
{
    delta_t = deltaT;
}

template <class Container>
const Container& YoshidaIntegrator<Container>::getContainer() const
{
    return *container;
}

template <class Container>
Container& YoshidaIntegrator<Container>::getContainer()
{
    return *container;
}

template <class Container>
void YoshidaIntegrator<Container>::setContainer(Container& container)
{
    YoshidaIntegrator<Container>::container = &container;
}

template <class Container>
const Acceleration::AccelerationAccumulator<Container>& YoshidaIntegrator<Container>::getAccumulator() const
{
    return *accumulator;
}

template <class Container>
Acceleration::AccelerationAccumulator<Container>& YoshidaIntegrator<Container>::getAccumulator()
{
    return *accumulator;
}

template <class Container>
void YoshidaIntegrator<Container>::setAccumulator(Acceleration::AccelerationAccumulator<Container>& accumulator)
{
    YoshidaIntegrator<Container>::accumulator = &accumulator;
}

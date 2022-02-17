//
// Created by Oliver on 13.05.21.
//

#include "Debris_test.h"
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/physics/AccelerationAccumulator.h"
#include "satellitePropagator/physics/Constants.h"
#include "satellitePropagator/physics/LeapFrogIntegrator.h"

/**
 * Initialize integrator, accumulator, and container with a single satellite in an equatorial orbit at given altitude.
 * @param satellite Satellite to be initialized. Only position and velocity are set.
 * @param altitude [km] Altitude of the satellite above ground.
 * @param perturbations Array of bools that specifies which perturbations are activated in the integrator.
 * @return tuple{orbitRadius, container, integrator, accumulator, keplerianOrbitDiscretizationSteps}
 */
auto initSystem(Debris::Debris& satellite, double altitude, const std::array<bool, 8>& perturbations)
{
    // calculate orbital parameters for satellite at 400km above earth
    const double orbitRadius = Physics::R_EARTH + altitude; // [km]
    const double keplerianOrbitSpeed = std::sqrt(Physics::GM_EARTH / orbitRadius); // [km/s]
    const double keplerianOrbitLength = 2. * M_PI * orbitRadius; // [km]
    const double keplerianOrbitAngularSpeed = keplerianOrbitSpeed / orbitRadius; // [1/s]
    const double keplerianOrbitPeriod = 2. * M_PI / keplerianOrbitAngularSpeed; // [s]
    // only move by at most one kilometer per discretization step
    const size_t keplerianOrbitDiscretizationSteps = std::ceil(keplerianOrbitLength);
    // choose deltaT so that end up at the point of origin
    const double deltaT = keplerianOrbitPeriod / static_cast<double>(keplerianOrbitDiscretizationSteps);

    auto container = std::make_unique<Debris::DebrisContainer<Debris::Debris>>();
    satellite.setPosition({ orbitRadius, 0., 0. });
    satellite.setVelocity({ 0., keplerianOrbitSpeed, 0. });
    container->addDebris(satellite);

    // create an integrator with all perturbations
    auto accumulator = std::make_unique<Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>>(
        perturbations,
        *container,
        0.0,
        nullptr);
    auto integrator = std::make_unique<LeapFrogIntegrator<Debris::DebrisContainer<Debris::Debris>>>(
        *container,
        *accumulator,
        deltaT);

    return std::make_tuple(orbitRadius, std::move(container), std::move(integrator), std::move(accumulator), keplerianOrbitDiscretizationSteps);
}

/**
 * Creates an evasivePreserving satellite and subject it to all perturbations.
 * 10 orbits are simulated and checked if its orbital height changes at the point of origin.
 */
TEST(DebrisTests, testEvasivePreserving)
{
    Debris::Debris satellite {};
    satellite.setActivityState(Debris::Debris::ActivityState::evasivePreserving);
    auto [orbitRadius, container, integrator, accumulator, keplerianOrbitDiscretizationSteps] = initSystem(
        satellite,
        400,
        std::array<bool, 8> {
            true, // KEP
            true, // J2
            true, // C22
            true, // S22
            true, // SOL
            true, // LUN
            true, // SRP
            true, // DRAG
        });

    // simulate 10 orbits
    const size_t numberOfOrbits = 10l;
    for (size_t i = 0; i < keplerianOrbitDiscretizationSteps * numberOfOrbits; ++i) {
        integrator->integrate(false);
        auto p = *container->begin();
        const auto& pPos = p.getPosition();
        ASSERT_EQ(pPos[2], 0.) << "Particle deviated from XY plane.";

        // calc sqrt because numbers grow large pretty quickly
        double altitude = std::sqrt(pPos[0] * pPos[0] + pPos[1] * pPos[1]);

        if (i % keplerianOrbitDiscretizationSteps == 0) {
            // check that height stays within 10cm each orbit
            ASSERT_NEAR(altitude, orbitRadius, 1e-4)
                << "Satellite deviated from keplerian orbit after "
                << (i / keplerianOrbitDiscretizationSteps)
                << " orbits";
        }
    }
}

/**
 * Creates an passive satellite and subject it only keplarian forces.
 * This should behave equivalent to the evasivePreserving particle.
 * 10 orbits are simulated and checked if its orbital height changes at the point of origin.
 */
TEST(DebrisTests, testPassive)
{
    Debris::Debris satellite {};
    satellite.setActivityState(Debris::Debris::ActivityState::passive);
    auto [orbitRadius, container, integrator, accumulator, keplerianOrbitDiscretizationSteps] = initSystem(
        satellite,
        400,
        std::array<bool, 8> {
            true, // KEP
            false, // J2
            false, // C22
            false, // S22
            false, // SOL
            false, // LUN
            false, // SRP
            false, // DRAG
        });

    // simulate 10 orbits
    const size_t numberOfOrbits = 10l;
    for (size_t i = 0; i < keplerianOrbitDiscretizationSteps * numberOfOrbits; ++i) {
        integrator->integrate(false);
        auto p = *container->begin();
        const auto& pPos = p.getPosition();
        ASSERT_EQ(pPos[2], 0.) << "Particle deviated from XY plane.";

        // calc sqrt because numbers grow large pretty quickly
        double altitude = std::sqrt(pPos[0] * pPos[0] + pPos[1] * pPos[1]);

        if (i % keplerianOrbitDiscretizationSteps == 0) {
            // check that height stays within 10cm each orbit
            ASSERT_NEAR(altitude, orbitRadius, 1e-4)
                << "Satellite deviated from keplerian orbit after "
                << (i / keplerianOrbitDiscretizationSteps)
                << " orbits";
        }
    }
}

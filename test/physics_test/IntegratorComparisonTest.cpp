/**
 * @file IntegratorComparisonTest.cpp
 * @author F. Gratl
 * @date 18.02.22
 */

#include "IntegratorComparisonTest.h"
#include "satellitePropagator/physics/AccelerationAccumulator.h"
#include "satellitePropagator/physics/YoshidaIntegrator.h"

TEST(IntegratorComparisonTest, YoshidaVsManually)
{
    Debris::DebrisContainer<Debris::Debris> container;

    Debris::Debris particle;
    particle.setPosition({ 8000.0, 0.0, 0.0 });
    particle.setVelocity({ 0.0, 8.0, 0.0 });
    container.addDebris(particle);

    const std::array<bool, 8> activatedForces { true, false, false, false, false, false, false, false };
    constexpr double deltaT { 10. };

    Acceleration::AccelerationAccumulator<decltype(container)> accelerationAccumulator(activatedForces, container, deltaT, nullptr);
    YoshidaIntegrator<decltype(container)> yoshidaIntegrator(container, accelerationAccumulator, deltaT);

    yoshidaIntegrator.integrate(false);

    const std::array<double, 3> expectedPos {8000, 80, 0}; // [km]
    const auto pos = container.begin()->getPosition();
    EXPECT_NEAR(pos[0], expectedPos[0], 0.5);
    EXPECT_NEAR(pos[1], expectedPos[1], 0.5);
    EXPECT_NEAR(pos[2], expectedPos[2], 0.);

    const std::array<double, 3> expectedVel {-0.06, 7.99, 0}; // [km/s]
    const auto vel = container.begin()->getVelocity();
    EXPECT_NEAR(vel[0], expectedVel[0], 0.01);
    EXPECT_NEAR(vel[1], expectedVel[1], 0.01);
    EXPECT_NEAR(vel[2], expectedVel[2], 0.);
}
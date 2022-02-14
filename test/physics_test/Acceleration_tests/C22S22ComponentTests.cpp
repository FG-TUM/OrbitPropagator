/**
 * @file C22S22ComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "C22S22ComponentTests.h"
#include "satellitePropagator/physics/AccelerationComponents/C22Component.h"
#include "satellitePropagator/physics/AccelerationComponents/C22S22Component.h"
#include "satellitePropagator/physics/AccelerationComponents/S22Component.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::C22S22Component::apply() is the same as adding up the results
 * of Acceleration::C22Component::apply() AND
 * Acceleration::S22Component::apply()
 */
TEST_F(C22S22ComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 12;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<std::array<double, 3>, num_debris> accelerations_3 {};
    std::array<double, 3> acc_total_dummy {};
    double t = 0;
    double c_term = std::cos(Physics::THETA_G + Physics::NU_EARTH * t);
    double s_term = std::sin(Physics::THETA_G + Physics::NU_EARTH * t);

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations_1[i] = Acceleration::C22Component::apply(container->getDebrisVector()[i], c_term,
            s_term);
        accelerations_2[i] = Acceleration::S22Component::apply(container->getDebrisVector()[i], c_term,
            s_term);
        accelerations_3[i] = Acceleration::C22S22Component::apply(container->getDebrisVector()[i], c_term,
            s_term);
    }

    // e-23 fails, but e-22 passes
    double abs_err = 1e-22;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations_1[i][0] + accelerations_2[i][0],
            accelerations_3[i][0], abs_err);
        EXPECT_NEAR(accelerations_1[i][1] + accelerations_2[i][1],
            accelerations_3[i][1], abs_err);
        EXPECT_NEAR(accelerations_1[i][2] + accelerations_2[i][2],
            accelerations_3[i][2], abs_err);
    }
}

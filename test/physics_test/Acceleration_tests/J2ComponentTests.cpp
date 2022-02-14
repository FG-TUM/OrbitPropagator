/**
 * @file J2ComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "J2ComponentTests.h"
#include "satellitePropagator/physics/AccelerationComponents/J2Component.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::J2Component::apply() is the same as using another
 * implementation
 */
TEST_F(J2ComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<double, 3> acc_total_dummy {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations_1[i] = Acceleration::J2Component::apply(container->getDebrisVector()[i]);
        calcJ2(container->getDebrisVector()[i], accelerations_2[i]);
    }

    // 10e-20 fails, but e-19 passes
    double abs_err = 1e-10;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations_1[i][0], accelerations_2[i][0], abs_err);
        EXPECT_NEAR(accelerations_1[i][1], accelerations_2[i][1], abs_err);
        EXPECT_NEAR(accelerations_1[i][2], accelerations_2[i][2], abs_err);
    }
}

/**
 * Tests if the acceleration calculated using
 * Acceleration::J2Component::apply() is the same as some hand calculated
 * values
 */
TEST_F(J2ComponentTests, EquilavelnceWIthPreCalculatedTest)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations {};
    std::array<double, 3> acc_total_dummy {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations[i] = Acceleration::J2Component::apply(container->getDebrisVector()[i]);
    }

    // e-21 fails, but e-20 passes
    double abs_err = 1e-20;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations[i][0], pre_calculated[i][0], abs_err);
        EXPECT_NEAR(accelerations[i][1], pre_calculated[i][1], abs_err);
        EXPECT_NEAR(accelerations[i][2], pre_calculated[i][2], abs_err);
    }
}

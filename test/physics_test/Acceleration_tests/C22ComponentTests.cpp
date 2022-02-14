/**
 * @file C22ComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "C22ComponentTests.h"
#include "satellitePropagator/physics/AccelerationComponents/C22Component.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::C22Component::apply() is the same as using another
 * implementation
 */
TEST_F(C22ComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 12;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<double, 3> acc_total_dummy {};
    double t = 0;
    double c_term = std::cos(Physics::THETA_G + Physics::NU_EARTH * t);
    double s_term = std::sin(Physics::THETA_G + Physics::NU_EARTH * t);

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations_1[i] = Acceleration::C22Component::apply(container->getDebrisVector()[i], c_term, s_term);
        calcC22(container->getDebrisVector()[i], t, accelerations_2[i]);
    }

    // e-23 fails, but e-22 passes
    double abs_err = 1e-22;
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
TEST_F(C22ComponentTests, EquilavelnceWIthPreCalculatedTest)
{
    GTEST_SKIP() << "Pre calculated values missing!";

    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations {};
    std::array<double, 3> acc_total_dummy {};
    double t = 0;
    double c_term = std::cos(Physics::THETA_G + Physics::NU_EARTH * t);
    double s_term = std::sin(Physics::THETA_G + Physics::NU_EARTH * t);

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations[i] = Acceleration::C22Component::apply(container->getDebrisVector()[i], c_term,
            s_term);
    }

    // 10e-22 fails, but e-21 passes
    double abs_err = 1e-0;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations[i][0], pre_calculated[i][0], abs_err);
        EXPECT_NEAR(accelerations[i][1], pre_calculated[i][1], abs_err);
        EXPECT_NEAR(accelerations[i][2], pre_calculated[i][2], abs_err);
    }
}

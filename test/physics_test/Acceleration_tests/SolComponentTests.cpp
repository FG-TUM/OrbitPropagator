/**
 * @file SolComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "SolComponentTests.h"
#include "satellitePropagator/physics/AccelerationComponents/SolComponent.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::SolComponent::apply() is the same as using a inefficient
 * implementation
 */
TEST_F(SolComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<double, 3> acc_total_dummy {};
    double t = 0.1;
    double d_ref;
    for (int j = 0; j < 10; ++j) {
        std::array<double, 6> sun_params = Acceleration::SolComponent::setUp(t * j);
        // calculate the acceleration for all particles using two different
        // functions
        for (int i = 0; i < num_debris; ++i) {
            accelerations_1[i] = Acceleration::SolComponent::apply(container->getDebrisVector()[i],
                d_ref,
                sun_params);
            calcSol(container->getDebrisVector()[i], t * j, accelerations_2[i]);
        }

        // 10e-20 fails, but e-19 passes
        double abs_err = 1e-20;
        for (int i = 0; i < num_debris; ++i) {
            EXPECT_NEAR(accelerations_1[i][0], accelerations_2[i][0], abs_err);
            EXPECT_NEAR(accelerations_1[i][1], accelerations_2[i][1], abs_err);
            EXPECT_NEAR(accelerations_1[i][2], accelerations_2[i][2], abs_err);
        }
    }
}

/**
 * Tests if the acceleration calculated using
 * Acceleration::SolComponent::apply() is the same as some and calculated
 * values
 */
TEST_F(SolComponentTests, EquilavelnceWIthPreCalculatedTest)
{
    GTEST_SKIP() << "Pre calculated values missing!";

//    const int num_debris = 9;
//    std::array<std::array<double, 3>, num_debris> accelerations {};
//    std::array<double, 3> acc_total_dummy {};
//
//    // calculate the acceleration for all particles using two different functions
//    for (int i = 0; i < num_debris; ++i) {
//        accelerations[i] = Acceleration::J2Component::apply(container->getDebrisVector()[i]);
//    }
//
//    // 10e-22 fails, but e-21 passes
//    double abs_err = 1e-0;
//    for (int i = 0; i < num_debris; ++i) {
//        EXPECT_NEAR(accelerations[i][0], pre_calculated[i][0], abs_err);
//        EXPECT_NEAR(accelerations[i][1], pre_calculated[i][1], abs_err);
//        EXPECT_NEAR(accelerations[i][2], pre_calculated[i][2], abs_err);
//    }
}

/**
 * Compares the result of calcSolParams() and
 * Acceleration::SolComponent::setUp()
 */
TEST_F(SolComponentTests, CompareSetupFunction)
{
    double t = 0.1;
    double abs_err = 1e-6;
    for (int i = 0; i < 10; ++i) {
        const std::array<double, 6> params_1 = Acceleration::SolComponent::setUp(t * i);
        const std::array<double, 6> params_2 = calcSolParams(t * i);
        ASSERT_NEAR(params_1[0], params_2[0], abs_err);
        ASSERT_NEAR(params_1[1], params_2[1], abs_err);
        ASSERT_NEAR(params_1[2], params_2[2], abs_err);
        ASSERT_NEAR(params_1[3], params_2[3], abs_err);
        ASSERT_NEAR(params_1[4], params_2[4], abs_err);
        ASSERT_NEAR(params_1[5], params_2[5], abs_err);
    }
}

/**
 * Compares the result of the Acceleration::SolComponent::apply() function
 * using the calcSolParams() function instead of
 * Acceleration::SolComponent::setUp()
 */
TEST_F(SolComponentTests, CompareAfterSetupCalculations)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<double, 3> acc_total_dummy {};
    double t = 0.1;
    double d_ref;
    for (int j = 0; j < 10; ++j) {
        std::array<double, 6> sun_params = calcSolParams(t * j);
        // calculate the acceleration for all particles using two different
        // functions
        for (int i = 0; i < num_debris; ++i) {
            accelerations_1[i] = Acceleration::SolComponent::apply(container->getDebrisVector()[i],
                d_ref,
                sun_params);
            calcSol(container->getDebrisVector()[i], t * j, accelerations_2[i]);
        }

        // 10e-20 fails, but e-19 passes
        double abs_err = 1e-20;
        for (int i = 0; i < num_debris; ++i) {
            EXPECT_NEAR(accelerations_1[i][0], accelerations_2[i][0], abs_err);
            EXPECT_NEAR(accelerations_1[i][1], accelerations_2[i][1], abs_err);
            EXPECT_NEAR(accelerations_1[i][2], accelerations_2[i][2], abs_err);
        }
    }
}

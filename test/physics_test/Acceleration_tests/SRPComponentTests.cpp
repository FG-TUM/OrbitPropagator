/**
 * @file SRPComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "SRPComponentTests.h"
#include "satellitePropagator/physics/AccelerationComponents/SRPComponent.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::SRPComponent::apply() is the same as using a inefficient
 * implementation
 */
TEST_F(SRPComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<double, 3> acc_total_dummy {};
    double t = 0.1;
    double d_ref = 0;
    for (int j = 0; j < 10; ++j) {
        std::array<double, 6> sun_params = Acceleration::SRPComponent::setUp(t * j);
        // calculate the acceleration for all particles using two different
        // functions
        for (int i = 0; i < container->getDebrisVector().size(); ++i) {
            accelerations_1[i] = Acceleration::SRPComponent::apply(container->getDebrisVector()[i],
                d_ref,
                sun_params);
            calcSRP(container->getDebrisVector()[i], t * j, accelerations_2[i]);
        }

        // result is identical
        double abs_err = 1e-16;
        for (int i = 0; i < container->getDebrisVector().size(); ++i) {
            EXPECT_NEAR(accelerations_1[i][0], accelerations_2[i][0], abs_err);
            EXPECT_NEAR(accelerations_1[i][1], accelerations_2[i][1], abs_err);
            EXPECT_NEAR(accelerations_1[i][2], accelerations_2[i][2], abs_err);
        }
    }
}

/**
 * Tests if the acceleration calculated using
 * Acceleration::SRPComponent::apply() is the same as some and calculated
 * values
 */
TEST_F(SRPComponentTests, EquilavelnceWIthPreCalculatedTest)
{
    GTEST_SKIP() << "Pre calculated values missing!";

//    const int num_debris = 9;
//    std::array<std::array<double, 3>, num_debris> accelerations {};
//    std::array<double, 3> acc_total_dummy {};
//
//    // calculate the acceleration for all particles using two different functions
//    for (int i = 0; i < container->getDebrisVector().size(); ++i) {
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

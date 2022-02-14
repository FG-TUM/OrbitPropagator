/**
 * @file KepComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "KepComponentTests.h"
#include "satellitePropagator/physics/AccelerationComponents/KepComponent.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::KepComponent::apply() is a function of
 * ||#Debris::Debris::position||
 */
TEST_F(KepComponentTests, RadialSymmetryTest)
{
    std::array<std::array<double, 3>, 8> accelerations {};
    for (auto& acc : accelerations) {
        acc[0] = 0;
        acc[1] = 0;
        acc[2] = 0;
    }
    // calculate the acceleration for all particles and accumulate the results in
    // accelerations[3]
    for (int i = 0; i < 3; ++i) {
        accelerations[i] = Acceleration::KepComponent::apply(container->getDebrisVector()[i]);
        accelerations[3][0] += accelerations[i][0];
        accelerations[3][1] += accelerations[i][1];
        accelerations[3][2] += accelerations[i][2];
    }
    for (int i = 3; i < 6; ++i) {
        accelerations[i + 1] = Acceleration::KepComponent::apply(container->getDebrisVector()[i]);
        accelerations[7][0] += accelerations[i + 1][0];
        accelerations[7][1] += accelerations[i + 1][1];
        accelerations[7][2] += accelerations[i + 1][2];
    }

    // only one axis non zero
    // only x should be a non zero value
    ASSERT_NE(accelerations[0][0], 0);
    ASSERT_EQ(accelerations[0][1], 0);
    ASSERT_EQ(accelerations[0][2], 0);
    // only y should be a non zero value
    ASSERT_EQ(accelerations[1][0], 0);
    ASSERT_NE(accelerations[1][1], 0);
    ASSERT_EQ(accelerations[1][2], 0);
    // only z should be a non zero value
    ASSERT_EQ(accelerations[2][0], 0);
    ASSERT_EQ(accelerations[2][1], 0);
    ASSERT_NE(accelerations[2][2], 0);
    // x of the first particle should be y of the second particle should be
    ASSERT_EQ(accelerations[0][0], accelerations[1][1]);
    // y of the first particle should be z of the second particle should be
    ASSERT_EQ(accelerations[1][1], accelerations[2][2]);
    // x, y and z of the array used as accumulator should all be the same and non
    // zero
    ASSERT_NE(accelerations[3][0], 0);
    ASSERT_EQ(accelerations[3][0], accelerations[3][1]);
    ASSERT_EQ(accelerations[3][1], accelerations[3][2]);

    // all axis non zero
    // x of the first particle should be y of the second particle should be
    ASSERT_EQ(accelerations[4][0], accelerations[5][1]);
    // y of the first particle should be z of the second particle should be
    ASSERT_EQ(accelerations[5][1], accelerations[6][2]);
    // x, y and z of the array used as accumulator should all be the same and non
    // zero
    ASSERT_NE(accelerations[7][0], 0);
    ASSERT_EQ(accelerations[7][0], accelerations[7][1]);
    ASSERT_EQ(accelerations[7][1], accelerations[7][2]);
}

/**
 * Tests if the acceleration calculated using
 * Acceleration::KepComponent::apply() is the same as using another
 * implementation
 */
TEST_F(KepComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 15;
    std::array<std::array<double, 3>, num_debris> accelerations_1 {};
    std::array<std::array<double, 3>, num_debris> accelerations_2 {};
    std::array<double, 3> acc_total_dummy {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations_1[i] = Acceleration::KepComponent::apply(container->getDebrisVector()[i]);
        calcKep(container->getDebrisVector()[i], accelerations_2[i]);
    }

    // e-17 fails, but e-16 passes
    double abs_err = 1e-16;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations_1[i][0], accelerations_2[i][0], abs_err);
        EXPECT_NEAR(accelerations_1[i][1], accelerations_2[i][1], abs_err);
        EXPECT_NEAR(accelerations_1[i][2], accelerations_2[i][2], abs_err);
    }
}

/**
 * Tests if the acceleration calculated using
 * Acceleration::KepComponent::apply() is the same as some hand calculated
 * values
 */
TEST_F(KepComponentTests, EquilavelnceWIthPreCalculatedTest)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations {};
    std::array<double, 3> acc_total_dummy {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 6; i < 6 + num_debris; ++i) {
        accelerations[i - 6] = Acceleration::KepComponent::apply(container->getDebrisVector()[i]);
    }

    // 1e-16 fails, but e-15 passes
    double abs_err = 1e-15;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations[i][0], pre_calculated[i][0], abs_err);
        EXPECT_NEAR(accelerations[i][1], pre_calculated[i][1], abs_err);
        EXPECT_NEAR(accelerations[i][2], pre_calculated[i][2], abs_err);
    }
}

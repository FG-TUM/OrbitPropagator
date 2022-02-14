/**
 * @file DragComponentTests.cpp
 * @author F. Gratl
 * @date 11.02.22
 */

#include "DragComponentTests.h"

/**
 * Tests if the acceleration calculated using
 * Acceleration::DragComponent::apply() is the same as using another
 * implementation
 */
TEST_F(DragComponentTests, CalculationEquivalenceTest)
{
    const int num_debris = 12;
    std::array<std::array<double, 3>, num_debris> dragAcc {};
    // expected acceleration
    std::array<std::array<double, 3>, num_debris> expected {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        dragAcc[i] = Acceleration::DragComponent::apply(container->getDebrisVector()[i]);
        calcDrag(container->getDebrisVector()[i], expected[i]);
    }

    // only error due to floating point arithmetic
    double abs_err = 5e-16;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(dragAcc[i][0], expected[i][0], abs_err);
        EXPECT_NEAR(dragAcc[i][1], expected[i][1], abs_err);
        EXPECT_NEAR(dragAcc[i][2], expected[i][2], abs_err);
    }
}

/**
 * Tests if the acceleration calculated using Acceleration::DragComponent::apply() is
 * in opposite direction of the velocity vector and within a reasonable magnitude.
 */
TEST_F(DragComponentTests, CheckDirectionAndMagnitude)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations[i] = Acceleration::DragComponent::apply(container->getDebrisVector()[i]);
    }

    for (int i = 0; i < num_debris; ++i) {
        // Check direction
        // aliases for readability
        const auto& initialVelocity = container->getDebrisVector()[i].getVelocity();
        const auto& calculatedDrag = accelerations[i];
        // expect cos to be near -1 => vectors point in almost opposite direction
        // drag doesn't push in exact opposing direction
        double abs_err = 0.002;
        EXPECT_NEAR(MathUtils::cosSimilarity(calculatedDrag, initialVelocity), -1., abs_err) << "Direction check failed for debris " << i;

        // Check magnitude
        const auto dragMagnitude = MathUtils::euclideanNorm(calculatedDrag);
        const auto velocityMagnitude = MathUtils::euclideanNorm(initialVelocity);
        EXPECT_GT(dragMagnitude, 0.) << "Drag is zero for debris " << i;
        // we expect drag to be significantly less than velocity
        // if input changes drastically this might be violated
        EXPECT_LE(dragMagnitude, velocityMagnitude * .03) << "Magnitude check failed for debris " << i;
    }
}

/**
 * Tests if the length magnitude of the acceleration calculated using
 * Acceleration::DragComponent::apply() is quadratic to the magnitude of the velocity vector
 */
TEST_F(DragComponentTests, CheckQuadraticToVelocity)
{
    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> v1_accelerations {};
    std::array<std::array<double, 3>, num_debris> v2_accelerations {};
    std::array<std::array<double, 3>, num_debris> v4_accelerations {};
    std::array<std::array<double, 3>, num_debris> v8_accelerations {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        Debris::Debris d = Debris::Debris(container->getDebrisVector()[i]);
        v1_accelerations[i] = Acceleration::DragComponent::apply(d);
        d.setVelocity({ d.getVelocity()[0] * 2., d.getVelocity()[1] * 2., d.getVelocity()[2] * 2. });
        v2_accelerations[i] = Acceleration::DragComponent::apply(d);
        d.setVelocity({ d.getVelocity()[0] * 2., d.getVelocity()[1] * 2., d.getVelocity()[2] * 2. });
        v4_accelerations[i] = Acceleration::DragComponent::apply(d);
        d.setVelocity({ d.getVelocity()[0] * 2., d.getVelocity()[1] * 2., d.getVelocity()[2] * 2. });
        v8_accelerations[i] = Acceleration::DragComponent::apply(d);
    }

    auto relError = [](auto actual, auto expected) {
        return std::abs((actual - expected) / expected);
    };

    // expect drag to rise roughly quadratically w.r.t. velocity (drag(2*v) == 4*drag(v))
    for (int i = 0; i < num_debris; ++i) {
        const auto v1norm = MathUtils::euclideanNorm(v1_accelerations[i]);
        const auto v2norm = MathUtils::euclideanNorm(v2_accelerations[i]);
        const auto v4norm = MathUtils::euclideanNorm(v4_accelerations[i]);
        const auto v8norm = MathUtils::euclideanNorm(v8_accelerations[i]);

        EXPECT_LE(relError(v1norm * 4, v2norm), 0.05) << "Failed for debris " << i;
        EXPECT_LE(relError(v2norm * 4, v4norm), 0.05) << "Failed for debris " << i;
        EXPECT_LE(relError(v4norm * 4, v8norm), 0.05) << "Failed for debris " << i;
    }
}

TEST_F(DragComponentTests, AtmosphericDesityTest)
{
    // assert the table is properly constructed
    auto prevAlt = 0.;
    for (const auto [alt, dens] : Acceleration::DragComponent::atmosphericDensityMap) {
        // altitude is strictly rising
        EXPECT_GT(alt, prevAlt);
        prevAlt = alt;
        // do not expect density to be falling monotonously

        // make sure that estimates match tables that are actually in the table
        const auto estimatedDensity = Acceleration::DragComponent::interpolateDensity(alt);
        EXPECT_NEAR(dens, estimatedDensity, 1e-15);
    }

    // Check that the linear interpolation produces values between interpolation points
    // Choose two points A and B, interpolate directly in the middle of it and compare to the midpoint
    const auto [altA, densA] = *std::next(Acceleration::DragComponent::atmosphericDensityMap.begin(), 2);
    const auto [altB, densB] = *std::next(Acceleration::DragComponent::atmosphericDensityMap.begin(), 3);
    const auto altMid = (altA + altB) / 2;
    const auto densMid = (densA + densB) / 2;
    const auto estimate = Acceleration::DragComponent::interpolateDensity(altMid);
    EXPECT_NEAR(estimate, densMid, 1e-15);
}

/**
 * Tests if the acceleration calculated using
 * Acceleration::DragComponent::apply() is the same as some hand calculated
 * values
 */
TEST_F(DragComponentTests, EquilavelnceWIthPreCalculatedTest)
{

    GTEST_SKIP() << "Pre calculated values missing!";

    const int num_debris = 9;
    std::array<std::array<double, 3>, num_debris> accelerations {};
    std::array<double, 3> acc_total_dummy {};

    // calculate the acceleration for all particles using two different functions
    for (int i = 0; i < num_debris; ++i) {
        accelerations[i] = Acceleration::DragComponent::apply(container->getDebrisVector()[i]);
    }

    double abs_err = 1e-20;
    for (int i = 0; i < num_debris; ++i) {
        EXPECT_NEAR(accelerations[i][0], pre_calculated[i][0], abs_err);
        EXPECT_NEAR(accelerations[i][1], pre_calculated[i][1], abs_err);
        EXPECT_NEAR(accelerations[i][2], pre_calculated[i][2], abs_err);
    }
}
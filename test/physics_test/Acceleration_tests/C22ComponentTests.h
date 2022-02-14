/**
 * @file C22ComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once
#include "gtest/gtest.h"
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/physics/Constants.h"

class C22ComponentTests : public ::testing::Test {
protected:
    std::shared_ptr<Debris::DebrisContainer<Debris::Debris>> container;
    std::array<std::array<double, 3>, 9> pre_calculated;

    virtual void SetUp()
    {
        container = std::make_shared<Debris::DebrisContainer<Debris::Debris>>();
        Debris::Debris d;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::array<double, 3> pos { 0, 0, 0 };
                pos[j] = (i + 2) * 3.5e3;
                d.setPosition(pos);
                container->addDebris(d);
            }
        }
        for (int i = 0; i < 3; ++i) {
            std::array<double, 3> pos { 0, 0, 0 };
            pos[i] = 5000;
            pos[(i + 1) % 3] = 4321;
            pos[(i + 2) % 3] = 3210;
            d.setPosition(pos);
            container->addDebris(d);
        }

        // calculated with wolfram alpha
        // TODO: calculate values
        pre_calculated[0] = { 0, 0, 0 };
        pre_calculated[1] = { 0, 0, 0 };
        pre_calculated[2] = { 0, 0, 0 };
        pre_calculated[3] = { 0, 0, 0 };
        pre_calculated[4] = { 0, 0, 0 };
        pre_calculated[5] = { 0, 0, 0 };
        pre_calculated[6] = { 0, 0, 0 };
        pre_calculated[7] = { 0, 0, 0 };
        pre_calculated[8] = { 0, 0, 0 };
    }

    template <class D>
    void calcC22(D& d, double t, std::array<double, 3>& acc_c22)
    {
        double gme = 3.986004407799724e+5;
        double re = 6378.1363;
        double c22 = 2.43914352398e-6;
        double theta = 280.4606;
        double nue = 4.178074622024230e-3;
        ASSERT_EQ(gme, Physics::GM_EARTH);
        ASSERT_EQ(re, Physics::R_EARTH);
        ASSERT_EQ(c22, Physics::C_22);
        ASSERT_EQ(theta, Physics::THETA_G);
        ASSERT_EQ(nue, Physics::NU_EARTH);
        double c_term = std::cos(theta + nue * t);
        double s_term = std::sin(theta + nue * t);
        double x = d.getPosition()[0] * c_term + d.getPosition()[1] * s_term;
        double y = -d.getPosition()[0] * s_term + d.getPosition()[1] * c_term;
        double z = d.getPosition()[2];
        double n1 = 5 * gme * re * re * std::sqrt(15) * c22 * (y * y - x * x);
        double n2 = gme * re * re * std::sqrt(15) * c22;
        double d1 = 2 * std::pow(x * x + y * y + z * z, 3.5);
        double d2 = std::pow(x * x + y * y + z * z, 2.5);
        double c22_x = (n1 * x) / d1 + (n2 * x) / d2;
        double c22_y = (n1 * y) / d1 - (n2 * y) / d2;
        double c22_z = (n1 * z) / d1;
        acc_c22[0] = c22_x * c_term - c22_y * s_term;
        acc_c22[1] = c22_x * s_term + c22_y * c_term;
        acc_c22[2] = c22_z;
    }
};

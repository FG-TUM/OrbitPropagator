/**
 * @file SRPComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once

#include <gtest/gtest.h>
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/physics/Constants.h"


class SRPComponentTests : public ::testing::Test {
protected:
    std::shared_ptr<Debris::DebrisContainer<Debris::Debris>> container;
    std::array<std::array<double, 3>, 9> pre_calculated;

    virtual void SetUp()
    {
        container = std::make_shared<Debris::DebrisContainer<Debris::Debris>>();
        Debris::Debris d;
        double t = 0;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::array<double, 3> pos { 0, 0, 0 };
                pos[j] = (i + 2) * 3.5e3;
                d.setPosition(pos);
                d.setAom(i * j + 0.1);
                container->addDebris(d);
            }
        }
        // calculated with wolfram alpha
        pre_calculated[0] = { 0, 0, 0 };
        pre_calculated[0] = { 0, 0, 0 };
        pre_calculated[2] = { 0, 0, 0 };
        pre_calculated[3] = { 0, 0, 0 };
        pre_calculated[4] = { 0, 0, 0 };
        pre_calculated[5] = { 0, 0, 0 };
        pre_calculated[6] = { 0, 0, 0 };
        pre_calculated[7] = { 0, 0, 0 };
        pre_calculated[8] = { 0, 0, 0 };
    }

    template <class D>
    void calcSRP(D& d, double t, std::array<double, 3>& acc_srp)
    {
        double p = 4.56e-3;
        double a = 1.49619e+8;
        double as = 149.619;
        double phis0 = 357.5256;
        double nus = 1.1407410259335311e-5;
        double o = 282.94;
        double e = 23.4392911;
        ASSERT_EQ(p, Physics::P_SRP);
        ASSERT_EQ(a, Physics::AU);
        ASSERT_EQ(as, Physics::AU_SCALED);
        ASSERT_EQ(phis0, Physics::PHI_SUN_0);
        ASSERT_EQ(nus, Physics::NU_SUN);
        ASSERT_EQ(o, Physics::OMEGA);
        ASSERT_EQ(e, Physics::EPSILON);

        double l = phis0 + nus * t;
        double r = as - 2.499 * std::cos(l * M_PIl / 180) - 0.021 * std::cos(2 * l * M_PIl / 180);
        double lambda = o + l + (6892.0 / 3600) * std::sin(l * M_PIl / 180) + (72.0 / 3600) * std::sin(2 * l * M_PIl / 180);
        double xs = r * std::cos(lambda * M_PIl / 180);
        double ys = r * std::sin(lambda * M_PIl / 180) * std::cos(e * M_PIl / 180);
        double zs = r * std::sin(lambda * M_PIl / 180) * std::sin(e * M_PIl / 180);
        xs *= 1e+6;
        ys *= 1e+6;
        zs *= 1e+6;

        double x = d.getPosition()[0];
        double y = d.getPosition()[1];
        double z = d.getPosition()[2];

        acc_srp[0] = (d.getAom() * p * a * a * (x - xs)) / std::pow((x - xs) * (x - xs) + (y - ys) * (y - ys) + (z - zs) * (z - zs), 1.5);
        acc_srp[1] = (d.getAom() * p * a * a * (y - ys)) / std::pow((x - xs) * (x - xs) + (y - ys) * (y - ys) + (z - zs) * (z - zs), 1.5);
        acc_srp[2] = (d.getAom() * p * a * a * (z - zs)) / std::pow((x - xs) * (x - xs) + (y - ys) * (y - ys) + (z - zs) * (z - zs), 1.5);
    }
};
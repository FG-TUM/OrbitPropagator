/**
 * @file SolComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once

#include <gtest/gtest.h>
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/physics/Constants.h"

class SolComponentTests : public ::testing::Test {
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
    void calcSol(D& d, double t, std::array<double, 3>& acc_sol)
    {
        double gms = 1.32712440018e+11;
        double phis0 = 357.5256;
        double nus = 1.1407410259335311e-5;
        double o = 282.94;
        double e = 23.4392911;
        ASSERT_EQ(gms, Physics::GM_SUN);
        ASSERT_EQ(phis0, Physics::PHI_SUN_0);
        ASSERT_EQ(nus, Physics::NU_SUN);
        ASSERT_EQ(o, Physics::OMEGA);
        ASSERT_EQ(e, Physics::EPSILON);

        double l = phis0 + nus * t;
        double r = 149.619 - 2.499 * std::cos(l * M_PIl / 180) - 0.021 * std::cos(2 * l * M_PIl / 180);
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
        double d1 = std::pow(
            (x - xs) * (x - xs) + (y - ys) * (y - ys) + (z - zs) * (z - zs), 1.5);
        double d2 = std::pow(xs * xs + ys * ys + zs * zs, 1.5);
        acc_sol[0] = -gms * ((x - xs) / d1 + xs / d2);
        acc_sol[1] = -gms * ((y - ys) / d1 + ys / d2);
        acc_sol[2] = -gms * ((z - zs) / d1 + zs / d2);
    }

    const std::array<double, 6> calcSolParams(double t)
    {
        std::array<double, 6> sol_params;
        double phis0 = 357.5256;
        double nus = 1.1407410259335311e-5;
        double o = 282.94;
        double e = 23.4392911;

        double l = phis0 + nus * t;
        double r = 149.619 - 2.499 * std::cos(l * M_PIl / 180) - 0.021 * std::cos(2 * l * M_PIl / 180);
        double lambda = o + l + (6892.0 / 3600) * std::sin(l * M_PIl / 180) + (72.0 / 3600) * std::sin(2 * l * M_PIl / 180);
        sol_params[0] = r * std::cos(lambda * M_PIl / 180);
        sol_params[1] = r * std::sin(lambda * M_PIl / 180) * std::cos(e * M_PIl / 180);
        sol_params[2] = r * std::sin(lambda * M_PIl / 180) * std::sin(e * M_PIl / 180);
        sol_params[0] *= 1e+6;
        sol_params[1] *= 1e+6;
        sol_params[2] *= 1e+6;
        double d2 = std::pow(sol_params[0] * sol_params[0] + sol_params[1] * sol_params[1] + sol_params[2] * sol_params[2],
            1.5);
        sol_params[3] = sol_params[0] / d2;
        sol_params[4] = sol_params[1] / d2;
        sol_params[5] = sol_params[2] / d2;
        return sol_params;
    }
};
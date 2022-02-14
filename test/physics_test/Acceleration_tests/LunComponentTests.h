/**
 * @file LunComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once
#include <gtest/gtest.h>
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/physics/Constants.h"

class LunComponentTests : public ::testing::Test {
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
        // calculated with wolfram alpha
        pre_calculated[0] = { 0, 0, 0 };
        pre_calculated[1] = { 0, 0, 0 };
        pre_calculated[2] = {
            0,
            0,

        };
        pre_calculated[3] = { 0, 0, 0 };
        pre_calculated[4] = { 0, 0, 0 };
        pre_calculated[5] = { 0, 0, 0 };
        pre_calculated[6] = { 0, 0, 0 };
        pre_calculated[7] = { 0, 0, 0 };
        pre_calculated[8] = { 0, 0, 0 };
    }

    const std::array<double, 6> calcMoonParams(double t)
    {
        const double e = 23.4392911;
        const double nus = 1.1407410259335311e-5;
        const double numa = 1.512151961904581e-4;
        const double nump = 1.2893925235125941e-6;
        const double nums = 6.128913003523574e-7;

        const double phi_m = nus * t;
        const double phi_m_a = numa * t;
        const double phi_m_p = nump * t;
        const double phi_m_s = nums * t;
        const double l_0 = phi_m_p + phi_m_a + 218.31617;
        const double l_m = phi_m_a + 134.96292;
        const double l1_m = phi_m + 357.52543;
        const double f_m = phi_m_p + phi_m_a + phi_m_s + 93.27283;
        const double d_m = phi_m_p + phi_m_a - phi_m + 297.85027;

        const double c_2dm = std::cos(2 * d_m * M_PIl / 180);
        const double c_lm = std::cos(l_m * M_PIl / 180);
        const double c_2lm = std::cos(2 * l_m * M_PIl / 180);
        const double c_l1m = std::cos(l1_m * M_PIl / 180);
        const double c_fm = std::cos(f_m * M_PIl / 180);
        const double c_2fm = std::cos(2 * f_m * M_PIl / 180);
        const double s_dm = std::sin(d_m * M_PIl / 180);
        const double s_2dm = std::sin(2 * d_m * M_PIl / 180);
        const double s_lm = std::sin(l_m * M_PIl / 180);
        const double s_2lm = std::sin(2 * l_m * M_PIl / 180);
        const double s_l1m = std::sin(l1_m * M_PIl / 180);
        const double s_fm = std::sin(f_m * M_PIl / 180);
        const double s_2fm = std::sin(2 * f_m * M_PIl / 180);

        double r_m = -152 * ((c_lm * c_l1m - s_lm * s_l1m) * c_2dm + (s_lm * c_l1m + c_lm * s_l1m) * s_2dm);
        r_m -= 171 * (c_lm * c_2dm - s_lm * s_2dm);
        r_m -= 205 * (c_l1m * c_2dm + s_l1m * s_2dm);
        r_m += 246 * (c_2lm * c_2dm + s_2lm * s_2dm);
        r_m -= 570 * c_2lm;
        r_m -= 2956 * c_2dm;
        r_m -= 3699 * (c_2dm * c_lm + s_2dm * s_lm);
        r_m -= 20905 * c_lm;
        r_m += 385000;

        double lambda_m = -(55.0 / 3600) * (s_2fm * c_2dm - c_2fm * s_2dm);
        lambda_m -= (110.0 / 3600) * (s_lm * c_l1m + c_lm * s_l1m);
        lambda_m -= (125.0 / 3600) * (s_dm);
        lambda_m += (148.0 / 3600) * (s_lm * c_l1m - c_lm * s_l1m);
        lambda_m -= (165.0 / 3600) * (s_l1m * c_2dm - c_l1m * s_2dm);
        lambda_m += (192.0 / 3600) * (s_lm * c_2dm + c_lm * s_2dm);
        lambda_m -= (206.0 / 3600) * ((s_lm * c_l1m + c_lm * s_l1m) * c_2dm - (c_lm * c_l1m - s_lm * s_l1m) * s_2dm);
        lambda_m -= (212.0 / 3600) * (s_2lm * c_2dm - c_2lm * s_2dm);
        lambda_m -= (412.0 / 3600) * (s_2fm);
        lambda_m -= (668.0 / 3600) * (s_l1m);
        lambda_m += (769.0 / 3600) * s_2lm;
        lambda_m += (2370.0 / 3600) * (s_2dm);
        lambda_m -= (4856.0 / 3600) * (s_lm * c_2dm - c_lm * s_2dm);
        lambda_m += (22640.0 / 3600) * s_lm;
        lambda_m += l_0;

        double beta_m = (11.0 / 3600) * ((s_fm * c_l1m - c_fm * s_l1m) * c_2dm - (c_fm * c_l1m + s_fm * s_l1m) * s_2dm);
        beta_m += (21.0 / 3600) * (s_fm * c_lm - c_fm * s_lm);
        beta_m -= (23.0 / 3600) * ((s_l1m * c_fm + c_l1m * s_fm) * c_2dm - (c_l1m * c_fm - s_l1m * s_fm) * s_2dm);
        beta_m -= (25.0 / 3600) * (s_fm * c_2lm - c_fm * s_2lm);
        beta_m -= (31.0 / 3600) * ((s_fm * c_lm - c_fm * s_lm) * c_2dm - (c_fm * c_lm + s_fm * s_lm) * s_2dm);
        beta_m += (44.0 / 3600) * ((s_lm * c_fm + c_lm * s_fm) * c_2dm - (c_lm * c_fm - s_lm * s_fm) * s_2dm);
        beta_m -= (526.0 / 3600) * (s_fm * c_2dm - c_fm * s_2dm);
        beta_m += (18520.0 / 3600) * std::sin((f_m + lambda_m - l_0 + (412.0 / 3600) * s_2fm + (541.0 / 3600) * s_l1m) * M_PIl / 180);

        std::array<double, 3> moon_pos;
        double c_term = std::cos(lambda_m * M_PIl / 180);
        double s_term = std::sin(lambda_m * M_PIl / 180);
        moon_pos[0] = c_term;
        moon_pos[1] = s_term;
        c_term = std::cos(beta_m * M_PIl / 180);
        s_term = std::sin(beta_m * M_PIl / 180);
        moon_pos[0] = moon_pos[0] * c_term;
        moon_pos[1] = moon_pos[1] * c_term;
        moon_pos[2] = s_term;
        c_term = std::cos(e * M_PIl / 180);
        s_term = std::sin(e * M_PIl / 180);

        std::array<double, 6> moon_params = {
            moon_pos[0],
            c_term * moon_pos[1] - s_term * moon_pos[2],
            s_term * moon_pos[1] + c_term * moon_pos[2],
            0,
            0,
            0
        };
        moon_params[0] = r_m * moon_params[0];
        moon_params[1] = r_m * moon_params[1];
        moon_params[2] = r_m * moon_params[2];
        double d2 = moon_params[0] * moon_params[0] + moon_params[1] * moon_params[1] + moon_params[2] * moon_params[2];
        d2 = 1 / std::sqrt(d2 * d2 * d2);
        moon_params[3] = moon_params[0] * d2;
        moon_params[4] = moon_params[1] * d2;
        moon_params[5] = moon_params[2] * d2;
        return moon_params;
    }

    template <class D>
    void calcLun(D& d, double t, std::array<double, 3>& acc_lun)
    {
        double gmm = 4.9028e+3;
        double e = 23.4392911;
        double nus = 1.1407410259335311e-5;
        double numa = 1.512151961904581e-4;
        double nump = 1.2893925235125941e-6;
        double nums = 6.128913003523574e-7;
        ASSERT_EQ(gmm, Physics::GM_MOON);
        ASSERT_EQ(e, Physics::EPSILON);
        ASSERT_EQ(nus, Physics::NU_SUN);
        ASSERT_EQ(numa, Physics::NU_MOON_A);
        ASSERT_EQ(nump, Physics::NU_MOON_P);
        ASSERT_EQ(nums, Physics::NU_MOON_S);

        double phi_m = nus * t;
        double phi_m_a = numa * t;
        double phi_m_p = nump * t;
        double phi_m_s = nums * t;
        double l_0 = phi_m_p + phi_m_a + 218.31617;
        double l_m = phi_m_a + 134.96292;
        double l1_m = phi_m + 357.52543;
        double f_m = phi_m_p + phi_m_a + phi_m_s + 93.27283;
        double d_m = phi_m_p + phi_m_a - phi_m + 297.85027;

        double r_m = -152 * std::cos((l_m + l1_m - 2 * d_m) * M_PIl / 180);
        r_m -= 171 * std::cos((l_m + 2 * d_m) * M_PIl / 180);
        r_m -= 205 * std::cos((l1_m - 2 * d_m) * M_PIl / 180);
        r_m += 246 * std::cos((2 * l_m - 2 * d_m) * M_PIl / 180);
        r_m -= 570 * std::cos((2 * l_m) * M_PIl / 180);
        r_m -= 2956 * std::cos((2 * d_m) * M_PIl / 180);
        r_m -= 3699 * std::cos((2 * d_m - l_m) * M_PIl / 180);
        r_m -= 20905 * std::cos(l_m * M_PIl / 180);
        r_m += 385000;

        double lambda_m = -(55.0 / 3600) * std::sin((2 * f_m - 2 * d_m) * M_PIl / 180);
        lambda_m -= (110.0 / 3600) * std::sin((l_m + l1_m) * M_PIl / 180);
        lambda_m -= (125.0 / 3600) * std::sin((d_m)*M_PIl / 180);
        lambda_m += (148.0 / 3600) * std::sin((l_m - l1_m) * M_PIl / 180);
        lambda_m -= (165.0 / 3600) * std::sin((l1_m - 2 * d_m) * M_PIl / 180);
        lambda_m += (192.0 / 3600) * std::sin((l_m + 2 * d_m) * M_PIl / 180);
        lambda_m -= (206.0 / 3600) * std::sin((l_m + l1_m - 2 * d_m) * M_PIl / 180);
        lambda_m -= (212.0 / 3600) * std::sin((2 * l_m - 2 * d_m) * M_PIl / 180);
        lambda_m -= (412.0 / 3600) * std::sin((2 * f_m) * M_PIl / 180);
        lambda_m -= (668.0 / 3600) * std::sin((l1_m)*M_PIl / 180);
        lambda_m += (769.0 / 3600) * std::sin((2 * l_m) * M_PIl / 180);
        lambda_m += (2370.0 / 3600) * std::sin((2 * d_m) * M_PIl / 180);
        lambda_m -= (4856.0 / 3600) * std::sin((l_m - 2 * d_m) * M_PIl / 180);
        lambda_m += (22640.0 / 3600) * std::sin((l_m)*M_PIl / 180);
        lambda_m += l_0;

        double beta_m = (11.0 / 3600) * std::sin((-l1_m + f_m - 2 * d_m) * M_PIl / 180);
        beta_m += (21.0 / 3600) * std::sin((-l_m + f_m) * M_PIl / 180);
        beta_m -= (23.0 / 3600) * std::sin((l1_m + f_m - 2 * d_m) * M_PIl / 180);
        beta_m -= (25.0 / 3600) * std::sin((-2 * l_m + f_m) * M_PIl / 180);
        beta_m -= (31.0 / 3600) * std::sin((-l_m + f_m - 2 * d_m) * M_PIl / 180);
        beta_m += (44.0 / 3600) * std::sin((l_m + f_m - 2 * d_m) * M_PIl / 180);
        beta_m -= (526.0 / 3600) * std::sin((f_m - 2 * d_m) * M_PIl / 180);
        beta_m += (18520.0 / 3600) * std::sin((f_m + lambda_m - l_0 + (412.0 / 3600) * std::sin((2 * f_m) * M_PIl / 180) + (541.0 / 3600) * std::sin((l1_m)*M_PIl / 180)) * M_PIl / 180);

        double s_lambda = std::sin(lambda_m * M_PIl / 180);
        double c_lambda = std::cos(lambda_m * M_PIl / 180);
        double s_beta = std::sin(beta_m * M_PIl / 180);
        double c_beta = std::cos(beta_m * M_PIl / 180);
        double s_e = std::sin(e * M_PIl / 180);
        double c_e = std::cos(e * M_PIl / 180);

        double xm = c_lambda * c_beta;
        double ym = c_e * s_lambda * c_beta - s_e * s_beta;
        double zm = c_e * s_beta + s_e * s_lambda * c_beta;

        xm *= r_m;
        ym *= r_m;
        zm *= r_m;

        double x = d.getPosition()[0];
        double y = d.getPosition()[1];
        double z = d.getPosition()[2];

        double d1 = std::pow(
            (x - xm) * (x - xm) + (y - ym) * (y - ym) + (z - zm) * (z - zm), 1.5);
        double d2 = std::pow(xm * xm + ym * ym + zm * zm, 1.5);

        acc_lun[0] = -gmm * (((x - xm) / d1) + (xm / d2));
        acc_lun[1] = -gmm * (((y - ym) / d1) + (ym / d2));
        acc_lun[2] = -gmm * (((z - zm) / d1) + (zm / d2));
    }
};
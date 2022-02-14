/**
 * @file J2ComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once

#include <gtest/gtest.h>
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"


class J2ComponentTests : public ::testing::Test {
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
        pre_calculated[0] = {
            -0.000010967387600710096587782776206071373077501857198479655576780, 0,
            0
        };
        pre_calculated[1] = {
            0, -0.000010967387600710096587782776206071373077501857198479655576780,
            0
        };
        pre_calculated[2] = {
            0, 0,
            0.0000219347752014201931755655524121427461550037143969593111535607
        };
        pre_calculated[3] = {
            -2.166397550757549943265733571569653941234934755255240607759e-6, 0, 0
        };
        pre_calculated[4] = {
            0, -2.166397550757549943265733571569653941234934755255240607759e-6, 0
        };
        pre_calculated[5] = {
            0, 0, 4.3327951015150998865314671431393078824698695105104812155181e-6
        };
        pre_calculated[6] = {
            -6.854617250443810367364235128794608173438660749049784735487e-7, 0, 0
        };
        pre_calculated[7] = {
            0, -6.854617250443810367364235128794608173438660749049784735487e-7, 0
        };
        pre_calculated[8] = {
            0, 0, 1.3709234500887620734728470257589216346877321498099569470975e-6
        };
    }

    template <class D>
    void calcJ2(D& d, std::array<double, 3>& acc_j2)
    {
        double gme = 3.986004407799724e+5;
        double re = 6378.1363;
        double c20 = -4.84165371736e-4;

        double x = d.getPosition()[0];
        double y = d.getPosition()[1];
        double z = d.getPosition()[2];
        double n1 = gme * re * re * std::sqrt(5) * c20;
        double d1 = 2 * std::sqrt(x * x + y * y + z * z);
        double d2 = (x * x + y * y + z * z) * (x * x + y * y + z * z);
        double n3 = 15 * (z * z);
        double d3 = (x * x + y * y + z * z) * (x * x + y * y + z * z) * (x * x + y * y + z * z);
        acc_j2[0] = ((n1 * x) / d1) * ((3 / d2) - (n3 / d3));
        acc_j2[1] = ((n1 * y) / d1) * ((3 / d2) - (n3 / d3));
        acc_j2[2] = ((n1 * z) / d1) * ((9 / d2) - (n3 / d3));
    }
};
/**
 * @file DragComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once

#include <gtest/gtest.h>
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"
#include "satellitePropagator/physics/Constants.h"
#include "satellitePropagator/physics/AccelerationComponents/DragComponent.h"

class DragComponentTests : public ::testing::Test {
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
                std::array<double, 3> vel { 4, 4, 4 };
                pos[j] = Physics::R_EARTH + (i + 2) * 100;
                d.setPosition(pos);
                d.setVelocity(vel);
                d.setBcInv(i * j + 0.1);
                container->addDebris(d);
            }
        }
        for (int i = 0; i < 3; ++i) {
            std::array<double, 3> pos { 0, 0, 0 };
            std::array<double, 3> vel { 4, 4, 4 };
            pos[i] = 5000;
            pos[(i + 1) % 3] = 4321;
            pos[(i + 2) % 3] = 3210;
            d.setPosition(pos);
            d.setVelocity(vel);
            d.setBcInv(i + 0.1);
            container->addDebris(d);
        }
    }
    template <class D>
    void calcDrag(D& d, std::array<double, 3>& acc_drag)
    {
        double re = 6378.1363;
        double p0 = 1.3e9;
        double h = 8.500;
        double oe = 7.292115e-5;
        ASSERT_EQ(re, Physics::R_EARTH);
        ASSERT_EQ(p0, Physics::P_GROUND);
        ASSERT_EQ(h, Physics::H_ATMOSPHERE);
        ASSERT_EQ(oe, Physics::ROT_ATMOSPHERE);

        double x = d.getPosition()[0];
        double y = d.getPosition()[1];
        double z = d.getPosition()[2];

        double pHeightOverGround = std::sqrt(x * x + y * y + z * z) - re;
        // old atmospheric model
        // double p = p0 * std::exp(-pHeightOverGround / h);
        // new table based atmospheric model
        double p = Acceleration::DragComponent::interpolateDensity(pHeightOverGround);

        double v_rel_x = d.getVelocity()[0] + oe * y;
        double v_rel_y = d.getVelocity()[0] - oe * x;
        double v_rel_z = d.getVelocity()[0];
        double v_rel_norm = MathUtils::euclideanNorm<std::array<double, 3>>({ v_rel_x, v_rel_y, v_rel_z });

        acc_drag[0] = -(p * v_rel_x * v_rel_norm * d.getBcInv()) / 2;
        acc_drag[1] = -(p * v_rel_y * v_rel_norm * d.getBcInv()) / 2;
        acc_drag[2] = -(p * v_rel_z * v_rel_norm * d.getBcInv()) / 2;
    }
};

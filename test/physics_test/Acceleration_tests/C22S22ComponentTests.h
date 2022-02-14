/**
 * @file C22S22ComponentTests.h
 * @author F. Gratl
 * @date 11.02.22
 */

#pragma once

#include <gtest/gtest.h>
#include "satellitePropagator/debris/Debris.h"
#include "satellitePropagator/debris/DebrisContainer.h"

class C22S22ComponentTests : public ::testing::Test {
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
    }
};
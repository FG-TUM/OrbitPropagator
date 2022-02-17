//
// Created by ubuntu on 18.05.21.
//

#pragma once
#include <iostream>
// project imports
#include <memory>

#include "satellitePropagator/debris/include.h"
#include "satellitePropagator/io/include.h"
#include "satellitePropagator/physics/include.h"


std::shared_ptr<Debris::DebrisContainer<Debris::Debris>> container;
std::shared_ptr<CommandLineInput> command_line;
std::shared_ptr<FileInput<Debris::DebrisContainer<Debris::Debris>>> file_input;
std::shared_ptr<FileOutput<Debris::DebrisContainer<Debris::Debris>>> file_output;
std::shared_ptr<Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>> accumulator;
std::shared_ptr<LeapFrogIntegrator<Debris::DebrisContainer<Debris::Debris>>> integrator;

int main(int argc, char** argv);

void initSimulation(int argc, char** argv);

void runSimulation();

void runThesisCalculations();
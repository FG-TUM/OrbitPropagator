//
// Created by Oliver on 12.05.21.
//

#include "debris_sim.h"

int main(int argc, char** argv)
{
    initSimulation(argc, argv);
    if (command_line->getRunThesisCalculations()){
        runThesisCalculations();
    }else{
        runSimulation();
    }
    return 0;
}


void initSimulation(int argc, char** argv)
{
    try {
        command_line = std::make_shared<CommandLineInput>(argc, argv);
        container = std::make_shared<Debris::DebrisContainer<Debris::Debris>>();
        file_input = std::make_shared<FileInput<Debris::DebrisContainer<Debris::Debris>>>(*container, command_line->getInputFilePath(),
            command_line->getInputFileType());
        file_output = std::make_shared<FileOutput<Debris::DebrisContainer<Debris::Debris>>>(*container, command_line->getOutputFilePath(),
            command_line->getOutputFileType(), file_input->getAccConfig());
        accumulator = std::make_shared<Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>>(
            file_input->getAccConfig(), *container, file_input->getStartT(), file_output.get());
        integrator = std::make_shared<LeapFrogIntegrator<Debris::DebrisContainer<Debris::Debris>>>(*container, *accumulator,
            file_input->getDeltaT());
    }
    // if something went wrong with the command line parsing
    catch (std::invalid_argument& e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

void runSimulation()
{
    double current_time = file_input->getStartT();
    if(!command_line->getRunThesisCalculations()){
        for (auto& d : container->getDebrisVector()) {
            std::cout << d.toString() << std::endl;
        }
    }
    int iteration = 0;
    double next_write_time = file_input->getStartT();
    while (current_time <= file_input->getEndT()) {
        iteration++;
        if (current_time >= next_write_time) {
            integrator->integrate(true);
            file_output->writeDebrisData(current_time);
            next_write_time += file_input->getWriteDeltaT();
        } else {
            integrator->integrate();
        }
        current_time += file_input->getDeltaT();
    }
    if(!command_line->getRunThesisCalculations()){
        for (auto& d : container->getDebrisVector()) {
            std::cout << d.toString() << std::endl;
        }
    }
}



void runThesisCalculations(){
    std::cout << "Starting thesis calculations" << std::endl;
    std::cout << "Run short delta t comparisons base 2" << std::endl;
    // read the delta t input file
    file_input->setInputFilePath("input/delta_t_short.txt");
    file_input->setInputFileType(InputFile::TXT);
    file_input->readDebrisData();
    // used to keep track of current delta t
    double delta_t = file_input->getDeltaT();
    // set accumulator values that won't change for the delta t simulations
    accumulator->setConfig(file_input->getAccConfig());
    for (int i = 1; delta_t >= std::pow(10,-9); ++i){
        std::cout << "Delta t: " << delta_t << std::endl << std::endl;
        // for each delta t we need new output files
        file_output = std::make_shared<FileOutput<Debris::DebrisContainer<Debris::Debris>>>(*container, "output/raw/delta_t_2_short_" + std::to_string(i) + ".csv",
                                                                                            OutputFile::CSV, file_input->getAccConfig());
        // set values for this delta t simulation
        container->cleanDebrisVector();
        file_input->readDebrisData();
        file_input->setDeltaT(delta_t);
        accumulator->setFileOutput(*file_output);
        accumulator->setT(file_input->getStartT());
        integrator->setDeltaT(delta_t);
        //run the simulation for the current delta t
        runSimulation();
        // update delta t to halt of it's value
        file_input->setDeltaT(delta_t/2);
        delta_t=file_input->getDeltaT();
    }

    std::cout << "Process csv data" << std::endl;
    delta_t = 1.;
    for (int i = 1; delta_t >= std::pow(10,-9); ++i){
        std::ifstream input("output/raw/delta_t_2_short_" + std::to_string(i) + ".csv");
        std::ofstream output("output/processed/delta_t_2_short_" + std::to_string(i) + ".csv");
        output << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        std::ifstream heyoka_input("output/raw/heyoka_delta_t_ld_sim.csv");
        std::array<long double, 3> vel{};
        std::array<long double, 3> heyoka_vel{};
        std::array<long double, 3> abs_err{};
        std::array<long double, 3> rel_err{};
        std::string val;
        if(input.is_open() && heyoka_input.is_open()){
            // skip column name line
            std::getline(input, val);
            std::getline(heyoka_input, val);
            // columns of processed
            output << "index,time,abs v err x,abs v err y,abs v err z,abs v err ave,rel v err x,rel v err y,rel v err z,rel v err ave" << std::endl;
            while(std::getline(heyoka_input,val,',') && std::getline(input,val,',')) {
                // copy index and time
                output << val << ",";
                std::getline(input, val, ',');
                output << val << ",";
                for (int i = 2; i < 6; ++i) {
                    // skip till column 7
                    std::getline(input, val, ',');
                }
                for (int i = 1; i < 42; ++i) {
                    // skip till column 43
                    std::getline(heyoka_input, val, ',');
                }

                //velocity
                std::getline(input, val, ',');
                std::istringstream(val) >> vel[0];
                std::getline(heyoka_input, val, ',');
                std::istringstream(val) >> heyoka_vel[0];

                std::getline(input, val, ',');
                std::istringstream(val) >> vel[1];
                std::getline(heyoka_input, val, ',');
                std::istringstream(val) >> heyoka_vel[1];

                std::getline(input, val, ',');
                std::istringstream(val) >> vel[2];
                std::getline(heyoka_input, val, ',');
                std::istringstream(val) >> heyoka_vel[2];

                abs_err[0] = std::abs(heyoka_vel[0]-vel[0]);
                abs_err[1] = std::abs(heyoka_vel[1]-vel[1]);
                abs_err[2] = std::abs(heyoka_vel[2]-vel[2]);
                IOUtils::to_ostream(abs_err,output,",",{"",","});
                output << (abs_err[0]+abs_err[1]+abs_err[2])/3. << ",";

                rel_err[0] = std::abs(abs_err[0]/heyoka_vel[0]);
                rel_err[1] = std::abs(abs_err[1]/heyoka_vel[1]);
                rel_err[2] = std::abs(abs_err[2]/heyoka_vel[2]);
                IOUtils::to_ostream(rel_err,output,",",{"",","});
                output << (rel_err[0]+rel_err[1]+rel_err[2])/3. << std::endl;

                // drop rest of the line
                std::getline(input,val);
                std::getline(heyoka_input,val);
            }
        }
        delta_t /= 2.;
    }

//    std::cout << "Run short delta t comparisons base 10" << std::endl;
//    // read the delta t input file
//    file_input->setInputFilePath("input/delta_t_short.txt");
//    file_input->setInputFileType(InputFile::TXT);
//    file_input->readDebrisData();
//    // used to keep track of current delta t
//    delta_t = file_input->getDeltaT();
//    // set accumulator values that won't change for the delta t simulations
//    accumulator->setConfig(file_input->getAccConfig());
//    for (int i = 1; delta_t >= std::pow(10,-8); ++i){
//        std::cout << "Delta t: " << delta_t << std::endl << std::endl;
//        // for each delta t we need new output files
//        file_output = std::make_shared<FileOutput<Debris::DebrisContainer<Debris::Debris>>>(*container, "output/raw/delta_t_10_short_" +
//        std::to_string(i) + ".csv",
//
//        OutputFile::CSV, file_input->getAccConfig());
//        // set values for this delta t simulation
//        container->cleanDebrisVector();
//        file_input->readDebrisData();
//        file_input->setDeltaT(delta_t);
//        accumulator->setFileOutput(*file_output);
//        accumulator->setT(file_input->getStartT());
//        integrator->setDeltaT(delta_t);
//        //run the simulation for the current delta t
//        runSimulation();
//        // update delta t to halt of it's value
//        file_input->setDeltaT(delta_t/10);
//        delta_t=file_input->getDeltaT();
//    }

//    std::cout << "Process csv data" << std::endl;
//    delta_t = 1.;
//    for (int i = 1; delta_t >= std::pow(10,-9); ++i){
//        std::ifstream input("output/raw/delta_t_10_short_" + std::to_string(i) + ".csv");
//        std::ofstream output("output/processed/delta_t_10_short_" + std::to_string(i) + ".csv");
//        output << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
//        std::ifstream heyoka_input("output/raw/heyoka_delta_t_ld_sim.csv");
//        std::array<long double, 3> vel{};
//        std::array<long double, 3> heyoka_vel{};
//        std::array<long double, 3> abs_err{};
//        std::array<long double, 3> rel_err{};
//        std::string val;
//        if(input.is_open() && heyoka_input.is_open()){
//            // skip column name line
//            std::getline(input, val);
//            std::getline(heyoka_input, val);
//            // columns of processed
//            output << "index,time,abs v err x,abs v err y,abs v err z,abs v err ave,rel v err x,rel v err y,rel v err z,rel v err ave" << std::endl;
//            while(std::getline(heyoka_input,val,',') && std::getline(input,val,',')) {
//
//                // copy index and time
//                output << val << ",";
//                std::getline(input, val, ',');
//                output << val << ",";
//                for (int i = 2; i < 6; ++i) {
//                    // skip till column 7
//                    std::getline(input, val, ',');
//                }
//                for (int i = 1; i < 42; ++i) {
//                    // skip till column 43
//                    std::getline(heyoka_input, val, ',');
//                }
//
//                //velocity
//                std::getline(input, val, ',');
//                std::istringstream(val) >> vel[0];
//                std::getline(heyoka_input, val, ',');
//                std::istringstream(val) >> heyoka_vel[0];
//
//                std::getline(input, val, ',');
//                std::istringstream(val) >> vel[1];
//                std::getline(heyoka_input, val, ',');
//                std::istringstream(val) >> heyoka_vel[1];
//
//                std::getline(input, val, ',');
//                std::istringstream(val) >> vel[2];
//                std::getline(heyoka_input, val, ',');
//                std::istringstream(val) >> heyoka_vel[2];
//
//                abs_err[0] = std::abs(heyoka_vel[0]-vel[0]);
//                abs_err[1] = std::abs(heyoka_vel[1]-vel[1]);
//                abs_err[2] = std::abs(heyoka_vel[2]-vel[2]);
//                IOUtils::to_ostream(abs_err,output,",",{"",","});
//                output << (abs_err[0]+abs_err[1]+abs_err[2])/3. << ",";
//
//                rel_err[0] = std::abs(abs_err[0]/heyoka_vel[0]);
//                rel_err[1] = std::abs(abs_err[1]/heyoka_vel[1]);
//                rel_err[2] = std::abs(abs_err[2]/heyoka_vel[2]);
//                IOUtils::to_ostream(rel_err,output,",",{"",","});
//                output << (rel_err[0]+rel_err[1]+rel_err[2])/3. << std::endl;
//
//                // drop rest of the line
//                std::getline(input,val);
//                std::getline(heyoka_input,val);
//            }
//        }
//        delta_t /= 10.;
//    }

//    std::cout << "Run long delta t comparisons base 10" << std::endl;
//    // read the delta t input file
//    file_input->setInputFilePath("input/delta_t.txt");
//    file_input->setInputFileType(InputFile::TXT);
//    file_input->readDebrisData();
//    // used to keep track of current delta t
//    delta_t = file_input->getDeltaT();
//    // set accumulator values that won't change for the delta t simulations
//    accumulator->setConfig(file_input->getAccConfig());
//    for (int i = 1; delta_t >= std::pow(10,-7); ++i){
//        std::cout << "Delta t: " << delta_t << std::endl << std::endl;
//        // for each delta t we need new output files
//        file_output = std::make_shared<FileOutput<Debris::DebrisContainer<Debris::Debris>>>(*container, "output/raw/delta_t_10_" + std::to_string(i) + ".csv",
//                                                                                            OutputFile::CSV, file_input->getAccConfig());
//        // set values for this delta t simulation
//        container->cleanDebrisVector();
//        file_input->readDebrisData();
//        file_input->setDeltaT(delta_t);
//        accumulator->setFileOutput(*file_output);
//        accumulator->setT(file_input->getStartT());
//        integrator->setDeltaT(delta_t);
//        //run the simulation for the current delta t
//        runSimulation();
//        // update delta t to halt of it's value
//        file_input->setDeltaT(delta_t/10);
//        delta_t=file_input->getDeltaT();
//    }

//    std::cout << "Process csv data" << std::endl;
//    delta_t = 1.;
//    for (int i = 1; delta_t >= std::pow(10,-9); ++i){
//        std::ifstream input("output/raw/delta_t_10_" + std::to_string(i) + ".csv");
//        std::ofstream output("output/processed/delta_t_10_" + std::to_string(i) + ".csv");
//        output << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
//        std::ifstream heyoka_input("output/raw/heyoka_delta_t_ld_sim.csv");
//        std::array<long double, 3> vel{};
//        std::array<long double, 3> heyoka_vel{};
//        std::array<long double, 3> abs_err{};
//        std::array<long double, 3> rel_err{};
//        std::string val;
//        if(input.is_open() && heyoka_input.is_open()){
//            // skip column name line
//            std::getline(input, val);
//            std::getline(heyoka_input, val);
//            // columns of processed
//            output << "index,time,abs v err x,abs v err y,abs v err z,abs v err ave,rel v err x,rel v err y,rel v err z,rel v err ave" << std::endl;
//            while(std::getline(heyoka_input,val,',') && std::getline(input,val,',')) {
//                // copy index and time
//                output << val << ",";
//                std::getline(input, val, ',');
//                output << val << ",";
//                for (int i = 2; i < 6; ++i) {
//                    // skip till column 7
//                    std::getline(input, val, ',');
//                }
//                for (int i = 1; i < 42; ++i) {
//                    // skip till column 43
//                    std::getline(heyoka_input, val, ',');
//                }
//
//                //velocity
//                std::getline(input, val, ',');
//                std::istringstream(val) >> vel[0];
//                std::getline(heyoka_input, val, ',');
//                std::istringstream(val) >> heyoka_vel[0];
//
//                std::getline(input, val, ',');
//                std::istringstream(val) >> vel[1];
//                std::getline(heyoka_input, val, ',');
//                std::istringstream(val) >> heyoka_vel[1];
//
//                std::getline(input, val, ',');
//                std::istringstream(val) >> vel[2];
//                std::getline(heyoka_input, val, ',');
//                std::istringstream(val) >> heyoka_vel[2];
//
//                abs_err[0] = std::abs(heyoka_vel[0]-vel[0]);
//                abs_err[1] = std::abs(heyoka_vel[1]-vel[1]);
//                abs_err[2] = std::abs(heyoka_vel[2]-vel[2]);
//                IOUtils::to_ostream(abs_err,output,",",{"",","});
//                output << (abs_err[0]+abs_err[1]+abs_err[2])/3. << ",";
//
//                rel_err[0] = std::abs(abs_err[0]/heyoka_vel[0]);
//                rel_err[1] = std::abs(abs_err[1]/heyoka_vel[1]);
//                rel_err[2] = std::abs(abs_err[2]/heyoka_vel[2]);
//                IOUtils::to_ostream(rel_err,output,",",{"",","});
//                output << (rel_err[0]+rel_err[1]+rel_err[2])/3. << std::endl;
//
//                // drop rest of the line
//                std::getline(input,val);
//                std::getline(heyoka_input,val);
//            }
//        }
//        delta_t /= 10.;
//    }
//    std::cout << "Calculations done" << std::endl;

//    file_input->setInputFilePath("input/orbit_comparison.txt");
//    file_input->setInputFileType(InputFile::TXT);
//    file_input->readDebrisData();
//    for (int day = 0; day < 366; ++day){
//        std::cout << "Delta t: " << delta_t << std::endl << std::endl;
//        // for each delta t we need new output files
//        file_output = std::make_shared<FileOutput<Debris::DebrisContainer<Debris::Debris>>>(*container, "output/raw/orbit_comparison" + std::to_string(day) + ".csv",
//                                                                                            OutputFile::CSV, file_input->getAccConfig());
//        // set values for this orbits simulation
//        container->cleanDebrisVector();
//        file_input->readDebrisData();
//        file_input->setStartT(day*24*60*60);
//        file_input->setEndT(file_input->getStartT() + file_input->getEndT());
//        accumulator->setFileOutput(*file_output);
//        accumulator->setT(file_input->getStartT());
//        integrator->setDeltaT(file_input->getDeltaT());
//        //run the simulation for the current delta t
//        runSimulation();
//    }

    std::ofstream out("output/raw/n_particles_runtime_multiple.csv");
    out << "n,runtime millis,millis per particle" << std::endl;
    file_input->setInputFilePath("input/runtime.txt");
    file_input->setInputFileType(InputFile::TXT);
    for(int j = 0;j < 10;++j) {
        container->cleanDebrisVector();
        file_input->readDebrisData();
        // run simulations with increasing number of particles
        for (int i = 0; i < 20; ++i) {
            accumulator->setT(file_input->getStartT());
            integrator->setDeltaT(file_input->getDeltaT());
            int n = container->getDebrisVector().size();
            std::cout << "n: " << n << std::endl;
            std::cout << "input start t: " << file_input->getStartT() << std::endl;
            std::cout << "accumulator t: " << accumulator->getT() << std::endl;
            std::cout << "input delta t: " << file_input->getDeltaT() << std::endl;
            std::cout << "integrator delta t: " << integrator->getDeltaT() << std::endl;
            std::cout << "input end t: " << file_input->getEndT() << std::endl;

            auto t1 = std::chrono::high_resolution_clock::now();
            runSimulation();
            auto t2 = std::chrono::high_resolution_clock::now();
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            std::cout << "Integrator runtime: " << ms_int.count() << " ms" << std::endl;
            out << n << "," << ms_int.count() << "," << ms_int.count() / (n * 1.) << std::endl;

            // double number of particles
            for (int j = 0; j < n; ++j) {
                container->addDebris(container->getDebrisVector()[j]);
            }
        }
    }
}

//
// Created by Oliver on 13.05.21.
//

#include "Integrator_test.h"

// setup taylor integrator and print it to check if the state etc. are correct
TEST_F(CompareWithHeyokaTests, showTaylorIntegrator)
{
    std::cout << *ta_total << std::endl;

    for(int i = 0; i < 8; ++i){
        std::cout << *ta_components[i] << std::endl;
    }
}

// compare calculated values of KepComponent
TEST_F(CompareWithHeyokaTests, compareKep) {
    // time step in seconds
    double delta_t = 1.0;
    // set some test debris values and add them to the Integrators
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        ds.push_back(d);
    }
    // loop over the debris data and compare calculations
    for (auto d : ds){
        // setup integrators
        i_components[Acceleration::KEP]->getDebris().cleanDebrisVector();
        i_components[Acceleration::KEP]->getDebris().addDebris(d);
        ta_components[Acceleration::KEP]->get_state_data()[0] = d.getPosition()[0];
        ta_components[Acceleration::KEP]->get_state_data()[1] = d.getPosition()[1];
        ta_components[Acceleration::KEP]->get_state_data()[2] = d.getPosition()[2];
        ta_components[Acceleration::KEP]->get_state_data()[3] = d.getVelocity()[0];
        ta_components[Acceleration::KEP]->get_state_data()[4] = d.getVelocity()[1];
        ta_components[Acceleration::KEP]->get_state_data()[5] = d.getVelocity()[2];
        // reset time values
        i_components[Acceleration::KEP]->setDeltaT(delta_t);
        i_components[Acceleration::KEP]->getAccumulator().setT(0.);
        ta_components[Acceleration::KEP]->set_time(0.);
        std::cout << *ta_components[Acceleration::KEP] << std::endl;
        // integrate time step
        i_components[Acceleration::KEP]->integrate();
        ta_components[Acceleration::KEP]->propagate_until(delta_t);
        std::cout << *ta_components[Acceleration::KEP] << std::endl;
        // compare result
        std::array<double,3> pos_i = i_components[Acceleration::KEP]->getDebris().getDebrisVector()[0].getPosition();
        std::array<double,3> vel_i = i_components[Acceleration::KEP]->getDebris().getDebrisVector()[0].getVelocity();
        std::array<double,3> pos_ta{ta_components[Acceleration::KEP]->get_state()[0],
                                    ta_components[Acceleration::KEP]->get_state()[1],
                                    ta_components[Acceleration::KEP]->get_state()[2]};
        std::array<double,3> vel_ta{ta_components[Acceleration::KEP]->get_state()[3],
                                    ta_components[Acceleration::KEP]->get_state()[4],
                                    ta_components[Acceleration::KEP]->get_state()[5]};
        IOUtils::to_ostream(pos_i, std::cout, ",", {"position integrator[","]\n"});
        IOUtils::to_ostream(pos_ta, std::cout, ",", {"position heyoka[","]\n"});
        IOUtils::to_ostream(vel_i, std::cout, ",", {"velocity integrator[","]\n"});
        IOUtils::to_ostream(vel_ta, std::cout, ",", {"velocity heyoka[","]\n"});
    }
}
// compare calculated values of J2Component
TEST_F(CompareWithHeyokaTests, compareJ2)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        ds.push_back(d);
    }

}
// compare calculated values of C22Component
TEST_F(CompareWithHeyokaTests, compareC22)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        ds.push_back(d);
    }

}
// compare calculated values of S22Component
TEST_F(CompareWithHeyokaTests, compareS22)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        ds.push_back(d);
    }

}
// compare calculated values of LunComponent
TEST_F(CompareWithHeyokaTests, compareLun)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        ds.push_back(d);
    }

}
// compare calculated values of SolComponent
TEST_F(CompareWithHeyokaTests, compareSol)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        ds.push_back(d);
    }

}
// compare calculated values of SRPComponent
TEST_F(CompareWithHeyokaTests, compareSRP)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        d.setAom(0.5);
        ds.push_back(d);
    }

}
// compare calculated values of DragComponent
TEST_F(CompareWithHeyokaTests, compareDrag)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        d.setBcInv(0.5);
        ds.push_back(d);
    }

}
// compare calculated values of all Components
TEST_F(CompareWithHeyokaTests, compareTotal)
{
    // set some test debris values
    std::vector<Debris::Debris> ds;
    Debris::Debris d;
    for (int i = 0; i < 3; ++i){
        d.setPosition({3500.*(i+2),0,0});
        d.setVelocity({10.,0,0});
        d.setAom(0.5);
        d.setBcInv(0.5);
        ds.push_back(d);
    }

}

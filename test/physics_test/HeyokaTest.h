//
// Created by Oliver on 13.05.21.
//

#pragma once
#ifdef HEYOKA_IS_PRESENT
#include "heyoka/heyoka.hpp"
#include "satellitePropagator/physics/Constants.h"
#include "satellitePropagator/physics/Integrator.h"
#include "gtest/gtest.h"

class CompareWithHeyokaTests : public ::testing::Test {
protected:
    //constants
    inline static constexpr double GMe = 3.986004407799724e+5;
    inline static constexpr double GMo = 1.32712440018e+11;
    inline static constexpr double GMm = 4.9028e+3;
    inline static constexpr double Re = 6378.1363;
    inline static constexpr double C20 = -4.84165371736e-4;
    inline static constexpr double C22 = 2.43914352398e-6;
    inline static constexpr double S22 = -1.40016683654e-6;
    inline static constexpr double theta_g = Physics::RAD_FACTOR * 280.4606;
    inline static constexpr double nu_e = Physics::RAD_FACTOR * (4.178074622024230e-3);
    inline static constexpr double nu_o = Physics::RAD_FACTOR * (1.1407410259335311e-5);
    inline static constexpr double nu_ma = Physics::RAD_FACTOR * (1.512151961904581e-4);
    inline static constexpr double nu_mp = Physics::RAD_FACTOR * (1.2893925235125941e-6);
    inline static constexpr double nu_ms = Physics::RAD_FACTOR * (6.128913003523574e-7);
    inline static constexpr double alpha_o = 1.49619e+8;
    inline static constexpr double epsilon = Physics::RAD_FACTOR * 23.4392911;
    inline static constexpr double phi_o = Physics::RAD_FACTOR * 357.5256;
    inline static constexpr double Omega_plus_w = Physics::RAD_FACTOR * 282.94;
    inline static constexpr double PSRP = 4.56e-3;
    inline static constexpr double p0 = 1.3;
    inline static constexpr double omega_e = 7.292115e-5;
    inline static constexpr double h = 8.5;

    // delta t
    inline static double delta_t = .1;
    inline static double start_t = 0.;
    inline static double end_t = 1000;

    //heyoka variables
    inline static heyoka::taylor_adaptive<double>* ta_total;
    inline static std::array<heyoka::taylor_adaptive<double>*, 8> ta_components;

    //heyoka integrator to extract single components
    inline static heyoka::taylor_adaptive<double>* ta_split;
    inline static heyoka::taylor_adaptive<long double>* ta_split_ld;

    // own integrators
    inline static Integrator<Debris::DebrisContainer<Debris::Debris>>* i_total;
    inline static Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>* aa_total;
    inline static std::array<Integrator<Debris::DebrisContainer<Debris::Debris>>*, 8> i_components;
    inline static std::array<Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>*, 8> aa_components;
    // output files
    inline static std::array<std::ofstream*, 8> ta_components_out;
    inline static std::array<FileOutput<Debris::DebrisContainer<Debris::Debris>>*, 8> i_components_out;
    inline static std::ofstream* ta_total_out;
    inline static FileOutput<Debris::DebrisContainer<Debris::Debris>>* i_total_out;
    inline static std::ofstream* ta_split_out;
    inline static std::ofstream* ta_split_ld_out;
    inline static std::array<int, 8> ta_components_line {};
    inline static int ta_total_line = 0;
    inline static int ta_split_line = 0;
    inline static int ta_split_ld_line = 0;
    inline static auto time_stamp = time(NULL);

    inline static void SetUpTestSuite()
    {
        //timestamp to create new data folder for each run
        std::filesystem::create_directory(std::filesystem::path(std::to_string(time_stamp)));
        //create heyoka variables
        std::array<heyoka::expression, 3> pos = heyoka::make_vars("X", "Y", "Z");
        std::array<heyoka::expression, 3> vel = heyoka::make_vars("VX", "VY", "VZ");

        //Sun's position
        auto lo = phi_o + nu_o * heyoka::time;
        auto lambda_o = Omega_plus_w + lo + Physics::RAD_FACTOR * ((6892. / 3600) * heyoka::sin(lo) + (72. / 3600) * heyoka::sin(2. * lo));
        auto ro = (149.619 - 2.499 * heyoka::cos(lo) - 0.021 * heyoka::cos(2. * lo)) * (1e6);

        auto Xo = ro * heyoka::cos(lambda_o);
        auto Yo = ro * heyoka::sin(lambda_o) * std::cos(epsilon);
        auto Zo = ro * heyoka::sin(lambda_o) * std::sin(epsilon);

        //Moon's position
        auto phi_m = nu_o * heyoka::time;
        auto phi_ma = nu_ma * heyoka::time;
        auto phi_mp = nu_mp * heyoka::time;
        auto phi_ms = nu_ms * heyoka::time;
        auto L0 = phi_mp + phi_ma + Physics::RAD_FACTOR * 218.31617;
        auto lm = phi_ma + Physics::RAD_FACTOR * 134.96292;
        auto llm = phi_m + Physics::RAD_FACTOR * 357.5256;
        auto Fm = phi_mp + phi_ma + phi_ms + Physics::RAD_FACTOR * 93.27283;
        auto Dm = phi_mp + phi_ma - phi_m + Physics::RAD_FACTOR * 297.85027;

        auto rm = 385000. - 20905. * heyoka::cos(lm) - 3699. * heyoka::cos(2. * Dm - lm) - 2956. * heyoka::cos(2. * Dm) - 570. * heyoka::cos(2. * lm) + 246. * heyoka::cos(2. * lm - 2. * Dm) - 205. * heyoka::cos(llm - 2. * Dm) - 171. * heyoka::cos(lm + 2. * Dm) - 152. * heyoka::cos(lm + llm - 2. * Dm);

        auto lambda_m = L0 + Physics::RAD_FACTOR * ((22640. / 3600) * heyoka::sin(lm) + (769. / 3600) * heyoka::sin(2. * lm) - (4856. / 3600) * heyoka::sin(lm - 2. * Dm) + (2370. / 3600) * heyoka::sin(2. * Dm) - (668. / 3600) * heyoka::sin(llm) - (412. / 3600) * heyoka::sin(2. * Fm) - (212. / 3600) * heyoka::sin(2. * lm - 2. * Dm) - (206. / 3600) * heyoka::sin(lm + llm - 2. * Dm) + (192. / 3600) * heyoka::sin(lm + 2. * Dm) - (165. / 3600) * heyoka::sin(llm - 2. * Dm) + (148. / 3600) * heyoka::sin(lm - llm) - (125. / 3600) * heyoka::sin(Dm) - (110. / 3600) * heyoka::sin(lm + llm) - (55. / 3600) * heyoka::sin(2. * Fm - 2. * Dm));

        auto Bm = Physics::RAD_FACTOR * ((18520. / 3600) * heyoka::sin(Fm + lambda_m - L0 + Physics::RAD_FACTOR * ((412. / 3600) * heyoka::sin(2. * Fm) + (541. / 3600) * heyoka::sin(llm))) - (526. / 3600) * heyoka::sin(Fm - 2. * Dm) + (44. / 3600) * heyoka::sin(lm + Fm - 2. * Dm) - (31. / 3600) * heyoka::sin(-lm + Fm - 2. * Dm) - (25. / 3600) * heyoka::sin(-2. * lm + Fm) - (23. / 3600) * heyoka::sin(llm + Fm - 2. * Dm) + (21. / 3600) * heyoka::sin(-lm + Fm) + (11. / 3600) * heyoka::sin(-llm + Fm - 2. * Dm));

        auto Xm = heyoka::cos(Bm) * heyoka::cos(lambda_m) * rm;
        auto Ym = -std::sin(epsilon) * heyoka::sin(Bm) * rm + std::cos(epsilon) * heyoka::cos(Bm) * heyoka::sin(lambda_m) * rm;
        auto Zm = std::cos(epsilon) * heyoka::sin(Bm) * rm + heyoka::cos(Bm) * std::sin(epsilon) * heyoka::sin(lambda_m) * rm;

        // distance to earth center
        auto magR2 = heyoka::pow(pos[0], 2.) + heyoka::pow(pos[1], 2.) + heyoka::pow(pos[2], 2.);

        //Earth's Keplerian terms
        auto fKepX = -GMe * pos[0] / (heyoka::pow(magR2, (3. / 2)));
        auto fKepY = -GMe * pos[1] / (heyoka::pow(magR2, (3. / 2)));
        auto fKepZ = -GMe * pos[2] / (heyoka::pow(magR2, (3. / 2)));

        //Earth's J2 terms
        auto J2term1 = GMe * (std::pow(Re, 2)) * std::sqrt(5) * C20 / (2. * heyoka::pow(magR2, (1. / 2)));
        auto J2term2 = 3. / (heyoka::pow(magR2, 2.));
        auto J2term3 = 15. * (heyoka::pow(pos[2], 2.)) / (heyoka::pow(magR2, 3.));
        auto fJ2X = J2term1 * pos[0] * (J2term2 - J2term3);
        auto fJ2Y = J2term1 * pos[1] * (J2term2 - J2term3);
        auto fJ2Z = J2term1 * pos[2] * (3. * J2term2 - J2term3);

        //Earth's C22 and S22 terms
        auto x = pos[0] * heyoka::cos(theta_g + nu_e * heyoka::time) + pos[1] * heyoka::sin(theta_g + nu_e * heyoka::time);
        auto y = -pos[0] * heyoka::sin(theta_g + nu_e * heyoka::time) + pos[1] * heyoka::cos(theta_g + nu_e * heyoka::time);
        auto z = pos[2];
        auto magr2 = heyoka::pow(x, 2.) + heyoka::pow(y, 2.) + heyoka::pow(z, 2.);

        auto C22term1 = 5 * GMe * (std::pow(Re, 2)) * std::sqrt(15) * C22 / (2. * heyoka::pow(magr2, (7. / 2)));
        auto C22term2 = GMe * (std::pow(Re, 2)) * std::sqrt(15) * C22 / (heyoka::pow(magr2, (5. / 2)));
        auto fC22x = C22term1 * x * (heyoka::pow(y, 2.) - heyoka::pow(x, 2.)) + C22term2 * x;
        auto fC22y = C22term1 * y * (heyoka::pow(y, 2.) - heyoka::pow(x, 2.)) - C22term2 * y;
        auto fC22z = C22term1 * z * (heyoka::pow(y, 2.) - heyoka::pow(x, 2.));

        auto S22term1 = 5 * GMe * (std::pow(Re, 2)) * std::sqrt(15) * S22 / (heyoka::pow(magr2, (7. / 2)));
        auto S22term2 = GMe * (std::pow(Re, 2)) * std::sqrt(15) * S22 / (heyoka::pow(magr2, (5. / 2)));
        auto fS22x = -S22term1 * (heyoka::pow(x, 2.)) * y + S22term2 * y;
        auto fS22y = -S22term1 * x * (heyoka::pow(y, 2.)) + S22term2 * x;
        auto fS22z = -S22term1 * x * y * z;

        auto fC22X = fC22x * heyoka::cos(theta_g + nu_e * heyoka::time) - fC22y * heyoka::sin(theta_g + nu_e * heyoka::time);
        auto fC22Y = fC22x * heyoka::sin(theta_g + nu_e * heyoka::time) + fC22y * heyoka::cos(theta_g + nu_e * heyoka::time);
        auto fC22Z = fC22z;

        auto fS22X = fS22x * heyoka::cos(theta_g + nu_e * heyoka::time) - fS22y * heyoka::sin(theta_g + nu_e * heyoka::time);
        auto fS22Y = fS22x * heyoka::sin(theta_g + nu_e * heyoka::time) + fS22y * heyoka::cos(theta_g + nu_e * heyoka::time);
        auto fS22Z = fS22z;

        //Sun's gravity
        auto magRo2 = heyoka::pow(Xo, 2.) + heyoka::pow(Yo, 2.) + heyoka::pow(Zo, 2.);
        auto magRRo2 = heyoka::pow((pos[0] - Xo), 2.) + heyoka::pow((pos[1] - Yo), 2.) + heyoka::pow((pos[2] - Zo), 2.);
        auto fSunX = -GMo * ((pos[0] - Xo) / (heyoka::pow(magRRo2, (3. / 2))) + Xo / (heyoka::pow(magRo2, (3. / 2))));
        auto fSunY = -GMo * ((pos[1] - Yo) / (heyoka::pow(magRRo2, (3. / 2))) + Yo / (heyoka::pow(magRo2, (3. / 2))));
        auto fSunZ = -GMo * ((pos[2] - Zo) / (heyoka::pow(magRRo2, (3. / 2))) + Zo / (heyoka::pow(magRo2, (3. / 2))));

        //Moon's gravity
        auto magRm2 = heyoka::pow(Xm, 2.) + heyoka::pow(Ym, 2.) + heyoka::pow(Zm, 2.);
        auto magRRm2 = heyoka::pow((pos[0] - Xm), 2.) + heyoka::pow((pos[1] - Ym), 2.) + heyoka::pow((pos[2] - Zm), 2.);
        auto fMoonX = -GMm * ((pos[0] - Xm) / (heyoka::pow(magRRm2, (3. / 2))) + Xm / (heyoka::pow(magRm2, (3. / 2))));
        auto fMoonY = -GMm * ((pos[1] - Ym) / (heyoka::pow(magRRm2, (3. / 2))) + Ym / (heyoka::pow(magRm2, (3. / 2))));
        auto fMoonZ = -GMm * ((pos[2] - Zm) / (heyoka::pow(magRRm2, (3. / 2))) + Zm / (heyoka::pow(magRm2, (3. / 2))));

        //Sun's radiation pressure (AOM is a heyoka parameter hy.par[0]. We
        //will be able to set it later without recompiling the integrator)
        auto SRPterm = heyoka::par[0] * PSRP * (std::pow(alpha_o, 2)) / (heyoka::pow(magRRo2, (3. / 2)));
        auto fSRPX = SRPterm * (pos[0] - Xo);
        auto fSRPY = SRPterm * (pos[1] - Yo);
        auto fSRPZ = SRPterm * (pos[2] - Zo);

        // atmospheric drag
        auto p = p0 * heyoka::exp(-(heyoka::sqrt(magR2) - Re) / h);
        auto v_rel_x = vel[0] + omega_e * pos[1];
        auto v_rel_y = vel[1] - omega_e * pos[0];
        auto v_rel_z = vel[2];
        // Drag due to the earths atmosphere (BCINV is a heyoka parameter hy.par[1])
        auto DragTerm = (heyoka::par[1] * p) / 2.;
        auto fDragX = -DragTerm * heyoka::pow(v_rel_x, 2.);
        auto fDragY = -DragTerm * heyoka::pow(v_rel_y, 2.);
        auto fDragZ = -DragTerm * heyoka::pow(v_rel_z, 2.);
        // EOM
        auto dXdt = vel[0];
        auto dYdt = vel[1];
        auto dZdt = vel[2];
        auto dVXdt = fKepX + fJ2X + fC22X + fS22X + fSunX + fMoonX + fSRPX + fDragX;
        auto dVYdt = fKepY + fJ2Y + fC22Y + fS22Y + fSunY + fMoonY + fSRPY + fDragY;
        auto dVZdt = fKepZ + fJ2Z + fC22Z + fS22Z + fSunZ + fMoonZ + fSRPZ + fDragZ;

        // set output accuracy
        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

        // Initial conditions (ic) -  can be changed later
        auto x0 = 0.;
        auto y0 = 0.;
        auto z0 = 0.;
        auto vx0 = 0.;
        auto vy0 = 0.;
        auto vz0 = 0.;

        auto t0 = 0.0;

        std::cout << "Compiling the Taylor integrator ... (this is done only once)" << std::endl;

        ta_total = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = dVXdt, heyoka::prime(vel[1]) = dVYdt, heyoka::prime(vel[2]) = dVZdt },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_total_out = new std::ofstream(std::to_string(time_stamp) + "/heyoka_total.csv");
        *ta_total_out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_total_out << "index,time,position x,position y,position z,position norm,velocity x, velocity y,velocity z,velocity norm" << std::endl;

        // calculates components separately

        std::array<heyoka::expression, 3> vel_kep = heyoka::make_vars("dXdtKepX", "dXdtKepY", "dXdtKepZ");
        auto vel_kep_x0 = 0.;
        auto vel_kep_y0 = 0.;
        auto vel_kep_z0 = 0.;

        std::array<heyoka::expression, 3> vel_j2 = heyoka::make_vars("dXdtJ2X", "dXdtJ2Y", "dXdtJ2Z");
        auto vel_j2_x0 = 0.;
        auto vel_j2_y0 = 0.;
        auto vel_j2_z0 = 0.;

        std::array<heyoka::expression, 3> vel_c22 = heyoka::make_vars("dXdtC22X", "dXdtC22Y", "dXdtC22Z");
        auto vel_c22_x0 = 0.;
        auto vel_c22_y0 = 0.;
        auto vel_c22_z0 = 0.;

        std::array<heyoka::expression, 3> vel_s22 = heyoka::make_vars("dXdtS22X", "dXdtS22Y", "dXdtS22Z");
        auto vel_s22_x0 = 0.;
        auto vel_s22_y0 = 0.;
        auto vel_s22_z0 = 0.;

        std::array<heyoka::expression, 3> vel_sun = heyoka::make_vars("dXdtSunX", "dXdtSunY", "dXdtSunZ");
        auto vel_sun_x0 = 0.;
        auto vel_sun_y0 = 0.;
        auto vel_sun_z0 = 0.;

        std::array<heyoka::expression, 3> vel_lun = heyoka::make_vars("dXdtLunX", "dXdtLunY", "dXdtLunZ");
        auto vel_lun_x0 = 0.;
        auto vel_lun_y0 = 0.;
        auto vel_lun_z0 = 0.;

        std::array<heyoka::expression, 3> vel_srp = heyoka::make_vars("dXdtSRPX", "dXdtSRPY", "dXdtSRPZ");
        auto vel_srp_x0 = 0.;
        auto vel_srp_y0 = 0.;
        auto vel_srp_z0 = 0.;

        std::array<heyoka::expression, 3> vel_drag = heyoka::make_vars("dXdtDragX", "dXdtDragY", "dXdtDragZ");
        auto vel_drag_x0 = 0.;
        auto vel_drag_y0 = 0.;
        auto vel_drag_z0 = 0.;

        auto dXXdt_split = vel[0] + vel_kep[0] + vel_j2[0] + vel_c22[0] + vel_s22[0] + vel_sun[0] + vel_lun[0] + vel_srp[0] + vel_drag[0];
        auto dXYdt_split = vel[1] + vel_kep[1] + vel_j2[1] + vel_c22[1] + vel_s22[1] + vel_sun[1] + vel_lun[1] + vel_srp[1] + vel_drag[1];
        auto dXZdt_split = vel[2] + vel_kep[2] + vel_j2[2] + vel_c22[2] + vel_s22[2] + vel_sun[2] + vel_lun[2] + vel_srp[2] + vel_drag[2];

        auto v_rel_x_split = dXXdt_split + omega_e * pos[1];
        auto v_rel_y_split = dXYdt_split - omega_e * pos[0];
        auto v_rel_z_split = dXZdt_split;
        auto fDragX_split = -DragTerm * heyoka::pow(v_rel_x_split, 2.);
        auto fDragY_split = -DragTerm * heyoka::pow(v_rel_y_split, 2.);
        auto fDragZ_split = -DragTerm * heyoka::pow(v_rel_z_split, 2.);

        ta_split = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXXdt_split, heyoka::prime(pos[1]) = dXYdt_split, heyoka::prime(pos[2]) = dXZdt_split,
                heyoka::prime(vel[0]) = heyoka::expression(0.), heyoka::prime(vel[1]) = heyoka::expression(0.), heyoka::prime(vel[2]) = heyoka::expression(0.),
                heyoka::prime(vel_kep[0]) = fKepX, heyoka::prime(vel_kep[1]) = fKepY, heyoka::prime(vel_kep[2]) = fKepZ,
                heyoka::prime(vel_j2[0]) = fJ2X, heyoka::prime(vel_j2[1]) = fJ2Y, heyoka::prime(vel_j2[2]) = fJ2Z,
                heyoka::prime(vel_c22[0]) = fC22X, heyoka::prime(vel_c22[1]) = fC22Y, heyoka::prime(vel_c22[2]) = fC22Z,
                heyoka::prime(vel_s22[0]) = fS22X, heyoka::prime(vel_s22[1]) = fS22Y, heyoka::prime(vel_s22[2]) = fS22Z,
                heyoka::prime(vel_sun[0]) = fSunX, heyoka::prime(vel_sun[1]) = fSunY, heyoka::prime(vel_sun[2]) = fSunZ,
                heyoka::prime(vel_lun[0]) = fMoonX, heyoka::prime(vel_lun[1]) = fMoonY, heyoka::prime(vel_lun[2]) = fMoonZ,
                heyoka::prime(vel_srp[0]) = fSRPX, heyoka::prime(vel_srp[1]) = fSRPY, heyoka::prime(vel_srp[2]) = fSRPZ,
                heyoka::prime(vel_drag[0]) = fDragX_split, heyoka::prime(vel_drag[1]) = fDragY_split, heyoka::prime(vel_drag[2]) = fDragZ_split },
            { x0, y0, z0,
                vx0, vy0, vz0,
                vel_kep_x0, vel_kep_y0, vel_kep_z0,
                vel_j2_x0, vel_j2_y0, vel_j2_z0,
                vel_c22_x0, vel_c22_y0, vel_c22_z0,
                vel_s22_x0, vel_s22_y0, vel_s22_z0,
                vel_sun_x0, vel_sun_y0, vel_sun_z0,
                vel_lun_x0, vel_lun_y0, vel_lun_z0,
                vel_srp_x0, vel_srp_y0, vel_srp_z0,
                vel_drag_x0, vel_drag_y0, vel_drag_z0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        ta_split_ld = new heyoka::taylor_adaptive<long double> {
            { heyoka::prime(pos[0]) = dXXdt_split, heyoka::prime(pos[1]) = dXYdt_split, heyoka::prime(pos[2]) = dXZdt_split,
                heyoka::prime(vel[0]) = heyoka::expression(0.), heyoka::prime(vel[1]) = heyoka::expression(0.), heyoka::prime(vel[2]) = heyoka::expression(0.),
                heyoka::prime(vel_kep[0]) = fKepX, heyoka::prime(vel_kep[1]) = fKepY, heyoka::prime(vel_kep[2]) = fKepZ,
                heyoka::prime(vel_j2[0]) = fJ2X, heyoka::prime(vel_j2[1]) = fJ2Y, heyoka::prime(vel_j2[2]) = fJ2Z,
                heyoka::prime(vel_c22[0]) = fC22X, heyoka::prime(vel_c22[1]) = fC22Y, heyoka::prime(vel_c22[2]) = fC22Z,
                heyoka::prime(vel_s22[0]) = fS22X, heyoka::prime(vel_s22[1]) = fS22Y, heyoka::prime(vel_s22[2]) = fS22Z,
                heyoka::prime(vel_sun[0]) = fSunX, heyoka::prime(vel_sun[1]) = fSunY, heyoka::prime(vel_sun[2]) = fSunZ,
                heyoka::prime(vel_lun[0]) = fMoonX, heyoka::prime(vel_lun[1]) = fMoonY, heyoka::prime(vel_lun[2]) = fMoonZ,
                heyoka::prime(vel_srp[0]) = fSRPX, heyoka::prime(vel_srp[1]) = fSRPY, heyoka::prime(vel_srp[2]) = fSRPZ,
                heyoka::prime(vel_drag[0]) = fDragX_split, heyoka::prime(vel_drag[1]) = fDragY_split, heyoka::prime(vel_drag[2]) = fDragZ_split },
            { x0, y0, z0,
                vx0, vy0, vz0,
                vel_kep_x0, vel_kep_y0, vel_kep_z0,
                vel_j2_x0, vel_j2_y0, vel_j2_z0,
                vel_c22_x0, vel_c22_y0, vel_c22_z0,
                vel_s22_x0, vel_s22_y0, vel_s22_z0,
                vel_sun_x0, vel_sun_y0, vel_sun_z0,
                vel_lun_x0, vel_lun_y0, vel_lun_z0,
                vel_srp_x0, vel_srp_y0, vel_srp_z0,
                vel_drag_x0, vel_drag_y0, vel_drag_z0 },
            heyoka::kw::time = t0,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_split_out = new std::ofstream(std::to_string(time_stamp) + "/heyoka_split.csv");
        ta_split_ld_out = new std::ofstream(std::to_string(time_stamp) + "/heyoka_split_ld.csv");
        *ta_split_out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_split_ld_out << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        *ta_split_out << "index,time,position x,position y,position z,position norm,vel_init x,vel_init y,vel_init z,vel_init norm,vel_kep x,vel_kep y,vel_kep z,vel_kep norm,vel_j2 x,vel_j2 y,vel_j2 z,vel_j2 norm,vel_c22 x,vel_c22 y,vel_c22 z,vel_c22 norm,vel_s22 x,vel_s22 y,vel_s22 z,vel_s22 norm,vel_sol x,vel_sol y,vel_sol z,vel_sol norm,vel_lun x,vel_lun y,vel_lun z,vel_lun norm,vel_srp x,vel_srp y,vel_srp z,vel_srp norm,vel_drag x,vel_drag y,vel_drag z,vel_drag norm,vel_total x,vel_total y,vel_total z,vel_total norm" << std::endl;
        *ta_split_ld_out << "index,time,position x,position y,position z,position norm,vel_init x,vel_init y,vel_init z,vel_init norm,vel_kep x,vel_kep y,vel_kep z,vel_kep norm,vel_j2 x,vel_j2 y,vel_j2 z,vel_j2 norm,vel_c22 x,vel_c22 y,vel_c22 z,vel_c22 norm,vel_s22 x,vel_s22 y,vel_s22 z,vel_s22 norm,vel_sol x,vel_sol y,vel_sol z,vel_sol norm,vel_lun x,vel_lun y,vel_lun z,vel_lun norm,vel_srp x,vel_srp y,vel_srp z,vel_srp norm,vel_drag x,vel_drag y,vel_drag z,vel_drag norm,vel_total x,vel_total y,vel_total z,vel_total norm" << std::endl;

        ta_components[Acceleration::KEP] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fKepX, heyoka::prime(vel[1]) = fKepY, heyoka::prime(vel[2]) = fKepZ },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::KEP] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_kep.csv");
        *ta_components_out[Acceleration::KEP] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::KEP] << "index,time,position x,position y,position z,position norm,vel_kep x,vel_kep y,vel_kep z,vel_kep norm" << std::endl;

        ta_components[Acceleration::J2] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fJ2X, heyoka::prime(vel[1]) = fJ2Y, heyoka::prime(vel[2]) = fJ2Z },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::J2] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_j2.csv");
        *ta_components_out[Acceleration::J2] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::J2] << "index,time,position x,position y,position z,position norm,vel_j2 x,vel_j2 y,vel_j2 z,vel_j2 norm" << std::endl;

        ta_components[Acceleration::C22] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fC22X, heyoka::prime(vel[1]) = fC22Y, heyoka::prime(vel[2]) = fC22Z },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::C22] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_c22.csv");
        *ta_components_out[Acceleration::C22] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::C22] << "index,time,position x,position y,position z,position norm,vel_c22 x,vel_c22 y,vel_c22 z,vel_c22 norm" << std::endl;

        ta_components[Acceleration::S22] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fS22X, heyoka::prime(vel[1]) = fS22Y, heyoka::prime(vel[2]) = fS22Z },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::S22] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_s22.csv");
        *ta_components_out[Acceleration::S22] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::S22] << "index,time,position x,position y,position z,position norm,vel_s22 x,vel_s22 y,vel_s22 z,vel_s22 norm" << std::endl;

        ta_components[Acceleration::SOL] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fSunX, heyoka::prime(vel[1]) = fSunY, heyoka::prime(vel[2]) = fSunZ },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::SOL] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_sol.csv");
        *ta_components_out[Acceleration::SOL] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::SOL] << "index,time,position x,position y,position z,position norm,vel_sol x,vel_sol y,vel_sol z,vel_sol norm" << std::endl;

        ta_components[Acceleration::LUN] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fMoonX, heyoka::prime(vel[1]) = fMoonY, heyoka::prime(vel[2]) = fMoonZ },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::LUN] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_lun.csv");
        *ta_components_out[Acceleration::LUN] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::LUN] << "index,time,position x,position y,position z,position norm,vel_lun x,vel_lun y,vel_lun z,vel_lun norm" << std::endl;

        ta_components[Acceleration::SRP] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fSRPX, heyoka::prime(vel[1]) = fSRPY, heyoka::prime(vel[2]) = fSRPZ },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::SRP] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_srp.csv");
        *ta_components_out[Acceleration::SRP] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::SRP] << "index,time,position x,position y,position z,position norm,vel_srp x,vel_srp y,vel_srp z,vel_srp norm" << std::endl;

        ta_components[Acceleration::DRAG] = new heyoka::taylor_adaptive<double> {
            { heyoka::prime(pos[0]) = dXdt, heyoka::prime(pos[1]) = dYdt, heyoka::prime(pos[2]) = dZdt, heyoka::prime(vel[0]) = fDragX, heyoka::prime(vel[1]) = fDragY, heyoka::prime(vel[2]) = fDragZ },
            { x0, y0, z0, vx0, vy0, vz0 },
            heyoka::kw::time = t0,
            heyoka::kw::tol = 1e-16,
            heyoka::kw::compact_mode = true
        };
        // prepare output file
        ta_components_out[Acceleration::DRAG] = new std::ofstream(std::to_string(time_stamp) + "/heyoka_drag.csv");
        *ta_components_out[Acceleration::DRAG] << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        *ta_components_out[Acceleration::DRAG] << "index,time,position x,position y,position z,position norm,vel_drag x,vel_drag y,vel_drag z,vel_drag norm" << std::endl;

        // setup own integrator for all components
        auto* container = new Debris::DebrisContainer<Debris::Debris>;
        std::array<bool, 8> config { true, true, true, true, true, true, true, true };

        i_total_out = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_total.csv"), OutputFile::Type::CSV, config);
        aa_total = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_total_out);
        i_total = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_total, delta_t);

        container = new Debris::DebrisContainer<Debris::Debris>;
        config = { false, false, false, false, false, false, false, false };
        config[Acceleration::KEP] = true;
        i_components_out[Acceleration::KEP] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_kep.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::KEP] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::KEP]);
        i_components[Acceleration::KEP] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::KEP], delta_t);
        config[Acceleration::KEP] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::J2] = true;
        i_components_out[Acceleration::J2] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_j2.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::J2] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::J2]);
        i_components[Acceleration::J2] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::J2], delta_t);
        config[Acceleration::J2] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::C22] = true;
        i_components_out[Acceleration::C22] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_c22.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::C22] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::C22]);
        i_components[Acceleration::C22] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::C22], delta_t);
        config[Acceleration::C22] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::S22] = true;
        i_components_out[Acceleration::S22] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_s22.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::S22] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::S22]);
        i_components[Acceleration::S22] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::S22], delta_t);
        config[Acceleration::S22] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::LUN] = true;
        i_components_out[Acceleration::LUN] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_lun.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::LUN] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::LUN]);
        i_components[Acceleration::LUN] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::LUN], delta_t);
        config[Acceleration::LUN] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::SOL] = true;
        i_components_out[Acceleration::SOL] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_sol.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::SOL] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::SOL]);
        i_components[Acceleration::SOL] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::SOL], delta_t);
        config[Acceleration::SOL] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::SRP] = true;
        i_components_out[Acceleration::SRP] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_srp.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::SRP] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::SRP]);
        i_components[Acceleration::SRP] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::SRP], delta_t);
        config[Acceleration::SRP] = false;

        container = new Debris::DebrisContainer<Debris::Debris>;
        config[Acceleration::DRAG] = true;
        i_components_out[Acceleration::DRAG] = new FileOutput<Debris::DebrisContainer<Debris::Debris>>(*container, std::filesystem::path(std::to_string(time_stamp) + "/integrator_drag.csv"), OutputFile::Type::CSV, config);
        aa_components[Acceleration::DRAG] = new Acceleration::AccelerationAccumulator<Debris::DebrisContainer<Debris::Debris>>(config, *container, start_t, i_components_out[Acceleration::DRAG]);
        i_components[Acceleration::DRAG] = new Integrator<Debris::DebrisContainer<Debris::Debris>>(*container, *aa_components[Acceleration::DRAG], delta_t);
        config[Acceleration::DRAG] = false;
    }

    template <class Container>
    void prepareRun(Integrator<Container>& i, heyoka::taylor_adaptive<double>& ta, typename Container::Particle_t& d)
    {
        i.getContainer().cleanDebrisVector();
        i.getContainer().addDebris(d);
        ta.get_state_data()[0] = d.getPosition()[0];
        ta.get_state_data()[1] = d.getPosition()[1];
        ta.get_state_data()[2] = d.getPosition()[2];
        ta.get_state_data()[3] = d.getVelocity()[0];
        ta.get_state_data()[4] = d.getVelocity()[1];
        ta.get_state_data()[5] = d.getVelocity()[2];
        // reset time values
        i.setDeltaT(delta_t);
        i.getAccumulator().setT(start_t);
        ta.set_time(start_t);
    }

    void showErrors(Integrator<Debris::DebrisContainer<Debris::Debris>>& i, heyoka::taylor_adaptive<double>& ta)
    {
        std::array<double, 3> pos_i = i.getContainer().getDebrisVector()[0].getPosition();
        std::array<double, 3> vel_i = i.getContainer().getDebrisVector()[0].getVelocity();
        std::array<double, 3> pos_ta { ta.get_state()[0],
            ta.get_state()[1],
            ta.get_state()[2] };
        std::array<double, 3> vel_ta { ta.get_state()[3],
            ta.get_state()[4],
            ta.get_state()[5] };
        IOUtils::to_ostream(pos_i, std::cout, ",", { "position integrator[", "]\n" });
        IOUtils::to_ostream(pos_ta, std::cout, ",", { "position heyoka[", "]\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::cosSimilarity(pos_ta, pos_i) }, std::cout, "",
            { "cosine similarity: ", "\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::euclideanDistance(pos_ta, pos_i) }, std::cout, "",
            { "euclidean distance: ", "\n" });
        IOUtils::to_ostream(MathUtils::absoluteError(pos_i, pos_ta), std::cout, ",", { "absolute error[", "]\n" });
        IOUtils::to_ostream(MathUtils::relativeError(pos_i, pos_ta), std::cout, ",", { "relative error[", "]\n" });
        IOUtils::to_ostream(vel_i, std::cout, ",", { "velocity integrator[", "]\n" });
        IOUtils::to_ostream(vel_ta, std::cout, ",", { "velocity heyoka[", "]\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::cosSimilarity(vel_ta, vel_i) }, std::cout, "",
            { "cosine similarity: ", "\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::euclideanDistance(vel_ta, vel_i) }, std::cout, "",
            { "euclidean distance: ", "\n" });
        IOUtils::to_ostream(MathUtils::absoluteError(vel_i, vel_ta), std::cout, ",", { "absolute error[", "]\n" });
        IOUtils::to_ostream(MathUtils::relativeError(vel_i, vel_ta), std::cout, ",", { "relative error[", "]\n" });
    }

    template <class D>
    void prepareRun(heyoka::taylor_adaptive<double>& split, heyoka::taylor_adaptive<double>& total, D& d)
    {
        split.get_state_data()[0] = d.getPosition()[0];
        split.get_state_data()[1] = d.getPosition()[1];
        split.get_state_data()[2] = d.getPosition()[2];
        split.get_state_data()[3] = d.getVelocity()[0];
        split.get_state_data()[4] = d.getVelocity()[1];
        split.get_state_data()[5] = d.getVelocity()[2];
        // we always start with 0 velocity
        for (int i = 6; i < split.get_state().size(); ++i) {
            split.get_state_data()[i] = 0;
        }
        total.get_state_data()[0] = d.getPosition()[0];
        total.get_state_data()[1] = d.getPosition()[1];
        total.get_state_data()[2] = d.getPosition()[2];
        total.get_state_data()[3] = d.getVelocity()[0];
        total.get_state_data()[4] = d.getVelocity()[1];
        total.get_state_data()[5] = d.getVelocity()[2];
        // reset time values
        split.set_time(start_t);
        total.set_time(start_t);
    }

    template <class D>
    void prepareLDRun(heyoka::taylor_adaptive<long double>& split, D& d)
    {
        split.get_state_data()[0] = d.getPosition()[0];
        split.get_state_data()[1] = d.getPosition()[1];
        split.get_state_data()[2] = d.getPosition()[2];
        split.get_state_data()[3] = d.getVelocity()[0];
        split.get_state_data()[4] = d.getVelocity()[1];
        split.get_state_data()[5] = d.getVelocity()[2];
        // we always start with 0 velocity
        for (int i = 6; i < split.get_state().size(); ++i) {
            split.get_state_data()[i] = 0;
        }
        // reset time values
        split.set_time(start_t);
    }

    void showErrors(heyoka::taylor_adaptive<double>& split, heyoka::taylor_adaptive<double>& total)
    {
        std::array<double, 3> pos_split { split.get_state()[0],
            split.get_state()[1],
            split.get_state()[2] };
        std::array<double, 3> vel_split {};
        // combine all components
        for (int i = 3; i < split.get_state().size(); i += 3) {
            vel_split[0] += split.get_state_data()[i];
            vel_split[1] += split.get_state_data()[i + 1];
            vel_split[2] += split.get_state_data()[i + 2];
        }
        std::array<double, 3> pos_ta { total.get_state()[0],
            total.get_state()[1],
            total.get_state()[2] };
        std::array<double, 3> vel_ta { total.get_state()[3],
            total.get_state()[4],
            total.get_state()[5] };
        IOUtils::to_ostream(pos_split, std::cout, ",", { "position integrator[", "]\n" });
        IOUtils::to_ostream(pos_ta, std::cout, ",", { "position heyoka[", "]\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::cosSimilarity(pos_ta, pos_split) }, std::cout, "",
            { "cosine similarity: ", "\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::euclideanDistance(pos_ta, pos_split) }, std::cout, "",
            { "euclidean distance: ", "\n" });
        IOUtils::to_ostream(MathUtils::absoluteError(pos_split, pos_ta), std::cout, ",", { "absolute error[", "]\n" });
        IOUtils::to_ostream(MathUtils::relativeError(pos_split, pos_ta), std::cout, ",", { "relative error[", "]\n" });
        IOUtils::to_ostream(vel_split, std::cout, ",", { "velocity integrator[", "]\n" });
        IOUtils::to_ostream(vel_ta, std::cout, ",", { "velocity heyoka[", "]\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::cosSimilarity(vel_ta, vel_split) }, std::cout, "",
            { "cosine similarity: ", "\n" });
        IOUtils::to_ostream(std::array<double, 1> { MathUtils::euclideanDistance(vel_ta, vel_split) }, std::cout, "",
            { "euclidean distance: ", "\n" });
        IOUtils::to_ostream(MathUtils::absoluteError(vel_split, vel_ta), std::cout, ",", { "absolute error[", "]\n" });
        IOUtils::to_ostream(MathUtils::relativeError(vel_split, vel_ta), std::cout, ",", { "relative error[", "]\n" });
    }

    void writeHeyokaState(Acceleration::AccelerationComponent type)
    {
        // one array that holds position and all velocity components in the order x,y,z,norm
        std::vector<double> output_vector;
        for (int i = 0; i < ta_components[type]->get_state().size() % 3; ++i) {
            // add x,y,z
            output_vector.push_back(ta_components[type]->get_state()[i]);
            output_vector.push_back(ta_components[type]->get_state()[i + 1]);
            output_vector.push_back(ta_components[type]->get_state()[i + 2]);
            output_vector.push_back(MathUtils::euclideanNorm(std::array<double, 3> { ta_components[type]->get_state()[i], ta_components[type]->get_state()[i + 1], ta_components[type]->get_state()[i + 2] }));
        }
        // output
        *ta_components_out[type] << ta_components_line[type]++ << ',';
        *ta_components_out[type] << ta_components[type]->get_time() << ',';
        IOUtils::to_ostream(output_vector, *ta_components_out[type]);
        *ta_components_out[type] << std::endl;
    }

    void writeSplitHeyokaState()
    {
        // one array that holds position and all velocity components in the order x,y,z,norm
        std::vector<double> output_vector;
        for (int i = 0; i < ta_split->get_state().size() / 3; ++i) {
            // add x,y,z
            output_vector.push_back(ta_split->get_state()[3 * i]);
            output_vector.push_back(ta_split->get_state()[3 * i + 1]);
            output_vector.push_back(ta_split->get_state()[3 * i + 2]);
            output_vector.push_back(MathUtils::euclideanNorm(std::array<double, 3> { ta_split->get_state()[3 * i], ta_split->get_state()[3 * i + 1], ta_split->get_state()[3 * i + 2] }));
        }
        std::array<double, 3> vel_split {};
        // combine all components
        for (int i = 3; i < ta_split->get_state().size(); i += 3) {
            vel_split[0] += ta_split->get_state_data()[i];
            vel_split[1] += ta_split->get_state_data()[i + 1];
            vel_split[2] += ta_split->get_state_data()[i + 2];
        }
        output_vector.push_back(vel_split[0]);
        output_vector.push_back(vel_split[1]);
        output_vector.push_back(vel_split[2]);
        output_vector.push_back(MathUtils::euclideanNorm(vel_split));

        // output
        *ta_split_out << ta_split_line++ << ',';
        *ta_split_out << ta_split->get_time() << ',';
        IOUtils::to_ostream(output_vector, *ta_split_out);
        *ta_split_out << std::endl;
    }

    void writeSplitLDHeyokaState()
    {
        // one array that holds position and all velocity components in the order x,y,z,norm
        std::vector<long double> output_vector;
        for (int i = 0; i < ta_split_ld->get_state().size() / 3; ++i) {
            // add x,y,z
            output_vector.push_back(ta_split_ld->get_state()[3 * i]);
            output_vector.push_back(ta_split_ld->get_state()[3 * i + 1]);
            output_vector.push_back(ta_split_ld->get_state()[3 * i + 2]);
            output_vector.push_back(MathUtils::euclideanNorm(std::array<long double, 3> { ta_split_ld->get_state()[3 * i], ta_split_ld->get_state()[3 * i + 1], ta_split_ld->get_state()[3 * i + 2] }));
        }
        std::array<double, 3> vel_split {};
        // combine all components
        for (int i = 3; i < ta_split_ld->get_state().size(); i += 3) {
            vel_split[0] += ta_split_ld->get_state_data()[i];
            vel_split[1] += ta_split_ld->get_state_data()[i + 1];
            vel_split[2] += ta_split_ld->get_state_data()[i + 2];
        }
        output_vector.push_back(vel_split[0]);
        output_vector.push_back(vel_split[1]);
        output_vector.push_back(vel_split[2]);
        output_vector.push_back(MathUtils::euclideanNorm(vel_split));

        // output
        *ta_split_ld_out << ta_split_ld_line++ << ',';
        *ta_split_ld_out << ta_split_ld->get_time() << ',';
        IOUtils::to_ostream(output_vector, *ta_split_ld_out);
        *ta_split_ld_out << std::endl;
    }
};
#endif

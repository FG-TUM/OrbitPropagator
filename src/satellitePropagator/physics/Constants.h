//
// Created by Oliver on 13.05.21.
//

#pragma once
#include <math.h>

// Define global constants
/**
 * @namespace Physics
 *
 * @brief Encapsulates global constants for calculations used in Acceleration
 * namespace
 */
namespace Physics {
/**
 * @brief Constant factor to get from degrees too radians
 */
inline constexpr double RAD_FACTOR = M_PIl / 180;
/**
 * @brief Gravitational constant * earth mass [km^3/s^2]
 *
 * Used in:
 * - Acceleration::KepComponent
 * - Acceleration::J2Component
 * - Acceleration::C22Component
 * - Acceleration::S22Component
 */
inline constexpr double GM_EARTH = 3.986004407799724e+5;
/**
 * @brief Gravitational constant * sun mass [km^3/s^2]
 *
 * Used in:
 * - Acceleration::SolComponent
 */
inline constexpr double GM_SUN = 1.32712440018e+11;
/**
 * @brief Gravitational constant * moon mass [km^3/s^2]
 *
 * Used in:
 * - Acceleration::LunComponent
 */
inline constexpr double GM_MOON = 4.9028e+3;
/**
 * @brief Radius of earth [km]
 *
 * Used in:
 * - Acceleration::J2Component
 * - Acceleration::C22Component
 * - Acceleration::S22Component
 */
inline constexpr double R_EARTH = 6378.1363;
/**
 * @brief Constant scalar factor
 *
 * Used in:
 * - Acceleration::J2Component
 */
inline constexpr double C_20 = -4.84165371736e-4;
/**
 * @brief Constant scalar factor
 *
 * Used in:
 * - Acceleration::C22Component
 */
inline constexpr double C_22 = 2.43914352398e-6;
/**
 * @brief Constant scalar factor
 *
 * Used in:
 * - Acceleration::S22Component
 */
inline constexpr double S_22 = -1.40016683654e-6;
/**
 * @brief Constant angle
 *
 * Used in:
 * - Acceleration::C22Component
 * - Acceleration::S22Component
 */
inline constexpr double THETA_G = 280.4606;
/**
 * @brief Angular velocity of earth rotation.
 *
 * Used in:
 * - Acceleration::C22Component
 * - Acceleration::S22Component
 */
inline constexpr double NU_EARTH = 4.178074622024230e-3;
/**
 * @brief Angular velocity of sun rotation
 *
 * Used in:
 * - Acceleration::SolComponent
 * - Acceleration::LunComponent
 */
inline constexpr double NU_SUN = 1.1407410259335311e-5;
/**
 * @brief Angular velocity of moon at point a
 *
 * Used in:
 * - Acceleration::LunComponent
 */
inline constexpr double NU_MOON_A = 1.512151961904581e-4;
/**
 * @brief Angular velocity of moon at point p
 *
 * Used in:
 * - Acceleration::LunComponent
 */
inline constexpr double NU_MOON_P = 1.2893925235125941e-6;
/**
 * @brief Angular velocity of moon at point s
 *
 * Used in:
 * - Acceleration::LunComponent
 */
inline constexpr double NU_MOON_S = 6.128913003523574e-7;
/**
 * @brief Astronomical unit
 *
 * Used in:
 * - Acceleration::SRPComponent
 */
inline constexpr double AU = 1.49619e+8;
/**
 * @brief Astronomical unit in 1e+6 km
 *
 * Used in:
 * - Acceleration::SolComponent
 */
inline constexpr double AU_SCALED = 149.619;
/**
 * @brief Constant angle
 *
 * Used in:
 * - Acceleration::SolComponent
 * - Acceleration::LunComponent
 */
inline constexpr double EPSILON = 23.4392911;
/**
 * @brief Constant angle
 *
 * Used in:
 * - Acceleration::SolComponent
 */
inline constexpr double PHI_SUN_0 = 357.5256;
/**
 * @brief Pressure of the radiation of the sun
 *
 * Used in:
 * - Acceleration::SRPComponent
 */
inline constexpr double P_SRP = 4.56e-3;
/**
 * @brief sum of longitude of Ascending node and argument of periapsis elements
 *
 * Used in:
 * - Acceleration::SolComponent
 */
inline constexpr double OMEGA = 282.94;
/**
 * @brief Atmospheric scale height [km].
 *
 * Used in:
 * - Acceleration::DragComponent
 */
inline constexpr double H_ATMOSPHERE = 8.5;
/**
 * @brief Rotation speed of earth around it's axis [rad/s]
 *
 * Used in:
 * - Acceleration::DragComponent
 */
inline constexpr double ROT_ATMOSPHERE = 7.292115e-5;
/**
 * @brief Atmospheric density at ground level [kg/km^3]
 *
 * Used in:
 * - Acceleration::DragComponent
 */
inline constexpr double P_GROUND = 1.3e9;
} // namespace Physics
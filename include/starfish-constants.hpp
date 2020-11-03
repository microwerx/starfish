#pragma once

#define _USE_MATH_DEFINES
#include <numeric>
#include <cmath>

//////////////////////////////////////////////////////////////////////
// M A T H E M A T I C A L   C O N S T A N T S ///////////////////////
//////////////////////////////////////////////////////////////////////

constexpr double SF_SQRT2 = 1.4142135623730950488016887242097;
constexpr double SF_SQRT3 = 1.7320508075688772935274463415059;
constexpr double SF_SQRT5 = 2.2360679774997896964091736687313;
constexpr double SF_ONEOVERSQRT2 = 0.70710678118654752440084436210485;
constexpr double SF_ONEOVERSQRT3 = 0.57735026918962576450914878050196;
constexpr double SF_ONEOVERSQRT5 = 0.44721359549995793928183473374626;
constexpr double SF_E = 2.7182818284590452353602874713527;
constexpr double SF_PHI = 1.618033988749894848204586834;
constexpr double SF_TWOPI = 6.283185307179586476925286766559;
constexpr double SF_PI = 3.1415926535897932384626433832795;
constexpr double SF_PIOVERTWO = 1.5707963267948966192313216916398;
constexpr double SF_PIOVERTHREE = 1.0471975511965977461542144610932;
constexpr double SF_PIOVERFOUR = 0.78539816339744830961566084581988;
constexpr double SF_RADIANS_TO_DEGREES = 57.295779513082320876798154814105;
constexpr double SF_EXTRA_RADIANS_TO_HOURS = 3.8197186342054880584532103209403;
constexpr double SF_EXTRA_RADIANS_TO_ARCMINS = 229.18311805232928350719261925642;
constexpr double SF_EXTRA_RADIANS_TO_ARCSECS = 13750.987083139757010431557155385;
constexpr double SF_DEGREES_TO_RADIANS = 0.01745329251994329576923690768489;
constexpr double SF_EXTRA_DEGREES_TO_HOURS = 0.06666666666666666666666666666667;
constexpr double SF_EXTRA_DEGREES_TO_ARCMINS = 4.0;
constexpr double SF_EXTRA_DEGREES_TO_ARCSECS = 240.0;
constexpr double SF_HOURS_TO_RADIANS = 0.26179938779914943653855361527329;
constexpr double SF_HOURS_TO_DEGREES = 15.0;
constexpr double SF_HOURS_TO_ARCMINS = 60.0;
constexpr double SF_HOURS_TO_ARCSECS = 3600.0;
constexpr double SF_ARCMINS_TO_RADIANS = 0.00436332312998582394230922692122;
constexpr double SF_ARCMINS_TO_DEGREES = 0.00416666666666666666666666666667;
constexpr double SF_ARCMINS_TO_HOURS = 0.01666666666666666666666666666667;
constexpr double SF_ARCMINS_TO_ARCSECS = 60.0;
constexpr double SF_ARCSECS_TO_DEGREES = 2.7777777777777777777777777777778e-4;
constexpr double SF_ARCSECS_TO_RADIANS = 4.8481368110953599358991410235795e-6;
constexpr double SF_ARCSECS_TO_HOURS = 2.7777777777777777777777777777778e-4;
constexpr double SF_ARCSECS_TO_ARCMINS = 0.01666666666666666666666666666667;
constexpr double SF_C = 299792458.0;            // speed of light
constexpr double SF_G = 6.6740831e-11;          // Newton's constant of gravitation
constexpr double SF_H = 6.62607004081e-34;      // Planck's h constant
constexpr double SF_HBAR = 1.05457180013e-34;   // Planck's h bar constant
constexpr double SF_EPSILON0 = 8.854187817e-12; // electric constant

constexpr double SF_LOG2E = 1.44269504088896340736;    // log2(e)
constexpr double SF_LOG10E = 0.434294481903251827651;  // log10(e)
constexpr double SF_LN2 = 0.693147180559945309417;     // ln(2)
constexpr double SF_LN10 = 2.30258509299404568402;     // ln(10)
constexpr double SF_PI_2 = 1.57079632679489661923;     // pi/2
constexpr double SF_PI_4 = 0.785398163397448309616;    // pi/4
constexpr double SF_1_PI = 0.318309886183790671538;    // 1/pi
constexpr double SF_2_PI = 0.636619772367581343076;    // 2/pi
constexpr double SF_1_2PI = 0.159154943;	           // 1/(2pi)
constexpr double SF_2_SQRTPI = 1.12837916709551257390; // 2/sqrt(pi)
constexpr double SF_SQRT1_2 = 0.707106781186547524401; // 1/sqrt(2)

constexpr float SF_F32_SQRT2 = 1.41421356f;
constexpr float SF_F32_SQRT3 = 1.73205080f;
constexpr float SF_F32_SQRT5 = 2.23606797f;
constexpr float SF_F32_ONEOVERSQRT2 = 0.70710678f;
constexpr float SF_F32_ONEOVERSQRT3 = 0.57735026f;
constexpr float SF_F32_ONEOVERSQRT5 = 0.44721359f;
constexpr float SF_F32_E = 2.71828182f;
constexpr float SF_F32_PHI = 1.61803398f;
constexpr float SF_F32_TWOPI = 6.28318530f;
constexpr float SF_F32_PI = 3.1415926535f;
constexpr float SF_F32_PIOVERTWO = 1.57079632f;
constexpr float SF_F32_PIOVERTHREE = 1.04719755f;
constexpr float SF_F32_PIOVERFOUR = 0.78539816f;
constexpr float SF_F32_RADIANS_TO_DEGREES = 57.29577951f;
constexpr float SF_F32_EXTRA_RADIANS_TO_HOURS = 3.81971863f;
constexpr float SF_F32_EXTRA_RADIANS_TO_ARCMINS = 229.18311805f;
constexpr float SF_F32_EXTRA_RADIANS_TO_ARCSECS = 13750.98708313f;
constexpr float SF_F32_DEGREES_TO_RADIANS = 0.01745329f;
constexpr float SF_F32_EXTRA_DEGREES_TO_HOURS = 0.06666666f;
constexpr float SF_F32_EXTRA_DEGREES_TO_ARCMINS = 4.0f;
constexpr float SF_F32_EXTRA_DEGREES_TO_ARCSECS = 240.0f;
constexpr float SF_F32_HOURS_TO_RADIANS = 0.26179938f;
constexpr float SF_F32_HOURS_TO_DEGREES = 15.0f;
constexpr float SF_F32_HOURS_TO_ARCMINS = 60.0f;
constexpr float SF_F32_HOURS_TO_ARCSECS = 3600.0f;
constexpr float SF_F32_ARCMINS_TO_RADIANS = 0.00436332f;
constexpr float SF_F32_ARCMINS_TO_DEGREES = 0.00416666f;
constexpr float SF_F32_ARCMINS_TO_HOURS = 0.01666666f;
constexpr float SF_F32_ARCMINS_TO_ARCSECS = 60.0f;
constexpr float SF_F32_ARCSECS_TO_DEGREES = 2.77777777e-4f;
constexpr float SF_F32_ARCSECS_TO_RADIANS = 4.84813681e-6f;
constexpr float SF_F32_ARCSECS_TO_HOURS = 2.77777777e-4f;
constexpr float SF_F32_ARCSECS_TO_ARCMINS = 0.01666666f;
constexpr float SF_F32_C = 299792458.0f;           // speed of light
constexpr float SF_F32_G = 6.6740831e-11f;         // Newton's constant of gravitation
constexpr float SF_F32_H = 6.62607004e-34f;        // Planck's h constant
constexpr float SF_F32_HBAR = 1.05457180e-34f;     // Planck's h bar constant
constexpr float SF_F32_EPSILON0 = 8.85418781e-12f; // electric constant

constexpr float SF_F32_LOG2E = 1.44269504f;    // log2(e)
constexpr float SF_F32_LOG10E = 0.43429448f;   // log10(e)
constexpr float SF_F32_LN2 = 0.69314718f;      // ln(2)
constexpr float SF_F32_LN10 = 2.30258509f;     // ln(10)
constexpr float SF_F32_PI_2 = 1.57079632f;     // pi/2
constexpr float SF_F32_PI_4 = 0.78539816f;     // pi/4
constexpr float SF_F32_1_PI = 0.31830988f;     // 1/pi
constexpr float SF_F32_2_PI = 0.63661977f;     // 2/pi
constexpr float SF_F32_1_2PI = 0.159154943f;   // 1/(2pi)
constexpr float SF_F32_2_SQRTPI = 1.12837916f; // 2/sqrt(pi)
constexpr float SF_F32_SQRT1_2 = 0.70710678f;  // 1/sqrt(2)

constexpr unsigned SF_INVALID_INDEX = 0xFFFFFFFF;

namespace Sf {
	template <typename T>
	constexpr T Radians(T angleInDegrees) noexcept {
		return (T)(angleInDegrees * SF_DEGREES_TO_RADIANS);
	}

	template <typename T>
	constexpr T Degrees(T angleInRadians) noexcept {
		return (T)(angleInRadians * SF_RADIANS_TO_DEGREES);
	}
}

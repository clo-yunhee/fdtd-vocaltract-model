#ifndef CONSTANTS_HH_
#define CONSTANTS_HH_

namespace vt {

namespace constants {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double sqrt_3 = 1.732050807568877293527446341505872366;

constexpr double rho = 1.149;  // Air density at 34°C [kg/m^3]
constexpr double c = 351.24;   // Speed of sound at 34°C [m/s]

}  // namespace constants

using namespace vt::constants;

}  // namespace vt

#endif  // CONSTANTS_HH_
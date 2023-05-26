#ifndef ROUTINES_BOUNDARY_INTERPOLATION_HH_
#define ROUTINES_BOUNDARY_INTERPOLATION_HH_

#include "../types.hh"

namespace vt {

array boundaryInterpolation(const array& tubeRadiusArray, double ds);

}

#endif  // ROUTINES_BOUNDARY_INTERPOLATION_HH_
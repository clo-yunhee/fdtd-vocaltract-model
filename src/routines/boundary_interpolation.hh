#ifndef ROUTINES_BOUNDARY_INTERPOLATION_HH_
#define ROUTINES_BOUNDARY_INTERPOLATION_HH_

#include "../types.hh"

namespace vt {

ArrayX<uint32_t> boundaryInterpolation(
    const Ref<ArrayX<uint32_t>> tubeRadiusArray, double ds);

}

#endif  // ROUTINES_BOUNDARY_INTERPOLATION_HH_
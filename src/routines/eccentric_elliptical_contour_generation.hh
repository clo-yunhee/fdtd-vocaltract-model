#ifndef ROUTINES_ECCENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
#define ROUTINES_ECCENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

array eccentricEllipticalContourGeneration(array& PV_N, uint32_t numSections,
                                           uint32_t     totalTubeLengthInCells,
                                           const array& tubeCummSectionLength,
                                           const array& ellipseAxisLenInfo,
                                           bool         boundaryInterpolation,
                                           StartInfo tubeStart, bool pmlSwitch,
                                           uint32_t pmlLayer, double ds);

}

#endif  // ROUTINES_ECCENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
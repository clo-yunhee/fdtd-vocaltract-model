#ifndef ROUTINES_CENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
#define ROUTINES_CENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

array centricEllipticalContourGeneration(array& PV_N, uint32_t numSections,
                                         uint32_t     totalTubeLengthInCells,
                                         const array& tubeCummSectionLength,
                                         const array& ellipseAxisLenInfo,
                                         bool         boundaryInterpolation,
                                         StartInfo tubeStart, double ds);

}

#endif  // ROUTINES_CENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
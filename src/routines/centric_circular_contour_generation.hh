#ifndef ROUTINES_CENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
#define ROUTINES_CENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

array centricCircularContourGeneration(array& PV_N, uint32_t numSections,
                                       uint32_t     totalTubeLengthInCells,
                                       const array& tubeSectionDiameterCells,
                                       const array& tubeCummSectionLength,
                                       bool         boundaryInterpolation,
                                       StartInfo tubeStart, double ds);

}

#endif  // ROUTINES_CENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
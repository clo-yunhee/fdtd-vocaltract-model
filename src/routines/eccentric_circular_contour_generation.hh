#ifndef ROUTINES_ECCENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
#define ROUTINES_ECCENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

array eccentricCircularContourGeneration(array& PV_N, uint32_t numSections,
                                         uint32_t     totalTubeLengthInCells,
                                         const array& tubeSectionDiameterCells,
                                         const array& tubeCummSectionLength,
                                         bool         boundaryInterpolation,
                                         StartInfo tubeStart, bool pmlSwitch,
                                         uint32_t pmlLayer, double ds);

}

#endif  // ROUTINES_ECCENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
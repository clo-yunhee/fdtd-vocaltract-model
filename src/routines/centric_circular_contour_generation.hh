#ifndef ROUTINES_CENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
#define ROUTINES_CENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

Array2X<uint32_t> centricCircularContourGeneration(
    Tensor4<double>& PV_N, uint32_t numSections,
    uint32_t                          totalTubeLengthInCells,
    const Ref<const ArrayX<uint32_t>> tubeSectionDiameterCells,
    const Ref<const ArrayXd> tubeCummSectionLength, bool boundaryInterpolation,
    StartInfo tubeStart, double ds);

}

#endif  // ROUTINES_CENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
#ifndef ROUTINES_ECCENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
#define ROUTINES_ECCENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

Array2X<uint32_t> eccentricCircularContourGeneration(
    Tensor4<double>& PV_N, uint32_t numSections,
    uint32_t                          totalTubeLengthInCells,
    const Ref<const ArrayX<uint32_t>> tubeSectionDiameterCells,
    const Ref<const ArrayXd> tubeCummSectionLength, bool boundaryInterpolation,
    StartInfo tubeStart, bool pmlSwitch, uint32_t pmlLayer, double ds);

}

#endif  // ROUTINES_ECCENTRIC_CIRCULAR_CONTOUR_GENERATION_HH_
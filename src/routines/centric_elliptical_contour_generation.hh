#ifndef ROUTINES_CENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
#define ROUTINES_CENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

Array2X<uint32_t> centricEllipticalContourGeneration(
    Tensor4<double>& PV_N, uint32_t numSections,
    uint32_t                           totalTubeLengthInCells,
    const Ref<const ArrayXd>           tubeCummSectionLength,
    const Ref<const Array2X<uint32_t>> ellipseAxisLenInfo,
    bool boundaryInterpolation, StartInfo tubeStart, double ds);

}

#endif  // ROUTINES_CENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
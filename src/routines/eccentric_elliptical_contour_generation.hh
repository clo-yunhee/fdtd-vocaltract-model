#ifndef ROUTINES_ECCENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
#define ROUTINES_ECCENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_

#include "../types.hh"

namespace vt {

Array2X<uint32_t> eccentricEllipticalContourGeneration(
    Tensor4<double>& PV_N, uint32_t numSections,
    uint32_t                           totalTubeLengthInCells,
    const Ref<const ArrayXd>           tubeCummSectionLength,
    const Ref<const Array2X<uint32_t>> ellipseAxisLenInfo,
    bool boundaryInterpolation, StartInfo tubeStart, bool pmlSwitch,
    uint32_t pmlLayer, double ds);

}

#endif  // ROUTINES_ECCENTRIC_ELLIPTICAL_CONTOUR_GENERATION_HH_
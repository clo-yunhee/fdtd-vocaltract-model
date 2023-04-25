#ifndef ROUTINES_FIND_CELL_TYPES_HH_
#define ROUTINES_FIND_CELL_TYPES_HH_

#include "../types.hh"

namespace vt {

ArrayXX<GridCellTypeInplane> findCellTypes(const Tensor4<double>& PV_N,
                                           uint32_t               tubeX);

}

#endif  // ROUTINES_FIND_CELL_TYPES_HH_
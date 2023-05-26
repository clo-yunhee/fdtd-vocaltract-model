#ifndef ROUTINES_TUBE_SEGMENT_CONNECTOR_HH_
#define ROUTINES_TUBE_SEGMENT_CONNECTOR_HH_

#include "../types.hh"

namespace vt {

void tubeSegmentConnector(array& PV_N, uint32_t prevGetRadius,
                          uint32_t currGetRadius, uint32_t tubeX);

}

#endif  // ROUTINES_TUBE_SEGMENT_CONNECTOR_HH_
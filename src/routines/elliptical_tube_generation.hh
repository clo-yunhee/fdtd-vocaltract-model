#ifndef ROUTINES_ELLIPTICAL_TUBE_GENERATION_HH_
#define ROUTINES_ELLIPTICAL_TUBE_GENERATION_HH_

#include "../types.hh"

namespace vt {

SimulationData ellipticalTubeGeneration(SimulationType simulationType,
                                        JunctionType junctionType, Vowel vowel,
                                        bool boundaryInterpolation,
                                        bool baffleSwitchFlag, bool pmlSwitch,
                                        uint32_t pmlLayer, MouthTermination rad,
                                        double ds);

}

#endif  // ROUTINES_ELLIPTICAL_TUBE_GENERATION_HH_
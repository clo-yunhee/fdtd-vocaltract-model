#ifndef PARAMS_HH_
#define PARAMS_HH_

#include "types.hh"

namespace vt {

struct Params {
    int              srateBase;  // Base sample rate (audio)
    uint32_t         srateMul;   // Sample rate multiplier (for simulation)
    double           dur;
    SimulationType   simulationType;
    Vowel            vowel;
    CrossSectionType crossSectionType;
    JunctionType     junctionType;
    SourceType       sourceType;
    double           LF_Rd;
    double           excitationF;
    double           srcAmplitude;
    bool             pmlSwitch;
    bool             baffleSwitch;
    MouthTermination rad;
    bool             boundaryInterpolation;
    bool             plotCells;
    bool             plotPressure;
    uint32_t         plotPressureStep;
    double           vis_Min;
    double           vis_Max;
    double           vis_Boundary;
};

}  // namespace vt

#endif  // PARAMS_HH_
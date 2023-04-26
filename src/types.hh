#ifndef TYPES_HH_
#define TYPES_HH_

#include <cstdint>

#include "eigen.hh"

namespace vt {

enum class SimulationType : uint32_t {
    OpenSpace,
    RegularTube,
    VowelSound,
};

enum class Vowel : uint32_t {
    None = 0,
    a,
    i,
    u,
    e,
    o,
    I,
};

enum class CrossSectionType : uint32_t {
    Circular,
    Elliptical,
    Square,
};

enum class JunctionType : uint32_t {
    Centric,
    Eccentric,
};

enum class SourceType : uint32_t {
    Sine,
    Gaussian,
    Impulse,
    GaussianNoise,
    VFModel,
    GFM_LF,
};

enum class MouthTermination : uint32_t {
    OpenMouth,
    DirichletBoundary,
};

enum GridCellType : uint32_t {
    cell_wall = 0,
    cell_air,
    cell_excitation,
    cell_pml0,
    cell_pml1,
    cell_pml2,
    cell_pml3,
    cell_pml4,
    cell_pml5,
    cell_dynamic,
    cell_dead,
    cell_noPressure,
    cell_head,
    // for book-keeping
    cell_numTypes,
};

inline bool operator==(const double& lhs, const GridCellType& rhs) {
    return (uint32_t)lhs == (uint32_t)rhs;
}

enum class GridCellTypeInplane : uint32_t {
    null = 0,
    inVTContour,
    outVTContour,
    onVTContour,
};

struct TubeProperties {
    ArrayXd sectionArea_incm2;
    double  totalLength;
    double  segmentLength;

    TubeProperties(std::initializer_list<double> p_sectionArea_incm2,
                   const double                  p_totalLength)
        : sectionArea_incm2(p_sectionArea_incm2.size() + 1),
          totalLength(p_totalLength),
          segmentLength(p_totalLength / p_sectionArea_incm2.size()) {
        std::copy(p_sectionArea_incm2.begin(), p_sectionArea_incm2.end(),
                  std::next(sectionArea_incm2.begin()));
    }
};

struct SimulationGridParams {
    uint32_t frameX;
    uint32_t frameY;
    uint32_t frameZ;
};

struct ListenerInfo {
    uint32_t listenerX;
    uint32_t listenerY;
    uint32_t listenerZ;
};

struct SourceInfo {
    uint32_t sourceX;
    uint32_t sourceY;
    uint32_t sourceZ;
};

struct StartInfo {
    uint32_t startX;
    uint32_t startY;
    uint32_t startZ;
};

struct EndInfo {
    uint32_t endX;
    uint32_t endY;
    uint32_t endZ;
};

struct SimulationData {
    SimulationGridParams gridParams;
    Tensor4<double>      PV_N;
    double               tubeStartArea;
    StartInfo            tubeStart;
    EndInfo              tubeEnd;
    uint32_t             totalTubeLengthInCells;
    Array2X<uint32_t>    currTubeSectionDiameterCells_SegmentCounter;
    ListenerInfo         listenerInfo;
    SourceInfo           sourceInfo;
};

}  // namespace vt

#endif  // TYPES_HH_
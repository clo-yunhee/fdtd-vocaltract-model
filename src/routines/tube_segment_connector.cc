#include "tube_segment_connector.hh"

#include "find_cell_types.hh"

void vt::tubeSegmentConnector(Tensor4<double>& PV_N, uint32_t prevGetRadius,
                              uint32_t currGetRadius, uint32_t tubeX) {
    // Store simulation grid size
    const auto& gridSize = dimensions(PV_N);

    // Move through the current grid planes and store the grid cells'
    // property in currGridPlaneProp
    const auto currGridPlaneProp = vt::findCellTypes(PV_N, tubeX);

    // Move through the previous grid planes and store the grid cells'
    // property in prevGridPlaneProp
    const auto prevGridPlaneProp = vt::findCellTypes(PV_N, tubeX - 1);

    // Based on the adjacent grid planes condition connect them
    if (currGetRadius > prevGetRadius) {
        for (uint32_t yCount = 1; yCount <= gridSize[1]; ++yCount) {
            for (uint32_t zCount = 1; zCount <= gridSize[3]; ++zCount) {
                if (currGridPlaneProp(yCount, zCount) ==
                        GridCellTypeInplane::inVTContour &&
                    prevGridPlaneProp(yCount, zCount) ==
                        GridCellTypeInplane::outVTContour) {
                    PV_N(yCount, tubeX - 1, zCount, 5) = vt::cell_wall;
                }
            }
        }
    } else {
        for (uint32_t yCount = 1; yCount <= gridSize[1]; ++yCount) {
            for (uint32_t zCount = 1; zCount <= gridSize[3]; ++zCount) {
                if (prevGridPlaneProp(yCount, zCount) ==
                        GridCellTypeInplane::inVTContour &&
                    currGridPlaneProp(yCount, zCount) ==
                        GridCellTypeInplane::outVTContour) {
                    PV_N(yCount, tubeX, zCount, 5) = vt::cell_wall;
                }
            }
        }
    }
}
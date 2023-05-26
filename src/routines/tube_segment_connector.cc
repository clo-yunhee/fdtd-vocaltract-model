#include "tube_segment_connector.hh"

#include "find_cell_types.hh"

void vt::tubeSegmentConnector(array& PV_N, uint32_t prevGetRadius,
                              uint32_t currGetRadius, uint32_t tubeX) {
    // Store simulation grid size
    const dim4 gridSize = PV_N.dims();

    // Move through the current grid planes and store the grid cells'
    // property in currGridPlaneProp
    const array currGridPlaneProp = vt::findCellTypes(PV_N, tubeX);

    // Move through the previous grid planes and store the grid cells'
    // property in prevGridPlaneProp
    const array prevGridPlaneProp = vt::findCellTypes(PV_N, tubeX - 1);

    // Based on the adjacent grid planes condition connect them
    if (currGetRadius > prevGetRadius) {
        for (uint32_t yCount = 0; yCount < gridSize[0]; ++yCount) {
            for (uint32_t zCount = 0; zCount < gridSize[2]; ++zCount) {
                if (allTrue<bool>(currGridPlaneProp(yCount, zCount) ==
                                  (float)GridCellTypeInplane::inVTContour) &&
                    allTrue<bool>(prevGridPlaneProp(yCount, zCount) ==
                                  (float)GridCellTypeInplane::outVTContour)) {
                    PV_N(yCount, tubeX - 1, zCount, 4) = vt::cell_wall;
                }
            }
        }
    } else {
        for (uint32_t yCount = 0; yCount < gridSize[0]; ++yCount) {
            for (uint32_t zCount = 0; zCount < gridSize[2]; ++zCount) {
                if (allTrue<bool>(prevGridPlaneProp(yCount, zCount) ==
                                  (float)GridCellTypeInplane::inVTContour) &&
                    allTrue<bool>(currGridPlaneProp(yCount, zCount) ==
                                  (float)GridCellTypeInplane::outVTContour)) {
                    PV_N(yCount, tubeX, zCount, 4) = vt::cell_wall;
                }
            }
        }
    }
}
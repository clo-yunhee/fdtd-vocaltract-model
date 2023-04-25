#include "circular_tube_generation.hh"

#include <cmath>

#include "../constants.hh"
#include "area_function.hh"
#include "centric_circular_contour_generation.hh"
#include "eccentric_circular_contour_generation.hh"
#include "find_cell_types.hh"

vt::SimulationData vt::circularTubeGeneration(
    SimulationType simulationType, JunctionType junctionType, Vowel vowel,
    bool boundaryInterpolation, bool baffleSwitchFlag, bool pmlSwitch,
    uint32_t pmlLayer, MouthTermination rad, double ds) {
    const double                 diameterMul = 1;
    ArrayXX<GridCellTypeInplane> gridPlaneProp;

    // Define tube geometry
    TubeProperties tube = vt::areaFunction(simulationType, junctionType, vowel);

    // Tube section area in m^2
    const auto tubeSectionArea_inm2 = tube.sectionArea_incm2 * 1e-2 * 1e-2;

    // Virtual mic/listener position near the mouth-end [inside the tube]
    const double   micXdistanceFromMouthEnd = 3e-3;  // 3 millimeters
    const uint32_t micXdistanceInCellsFromMouthEnd =
        (uint32_t)std::round(micXdistanceFromMouthEnd / ds);

    // Number of tube segments
    const uint32_t numSections = length(tubeSectionArea_inm2);

    // Calculate tube section diameters
    const auto tubeSectionOriginalDiameter =
        2 * sqrt(tubeSectionArea_inm2 / pi) * diameterMul;
    ArrayX<uint32_t> tubeSectionDiameterCells =
        round(tubeSectionOriginalDiameter / ds).cast<uint32_t>();

    // Cross-sectional area at the tube start end
    // I didn't use "tubeSectionArea_inm2" directly because I wanted to use
    // the derived area by using the derived diameter.
    const double tubeStartArea =
        pi * pow(tubeSectionOriginalDiameter(1) / 2, 2);

    // Change the tubeSectionDiameterCells to 1 if it contains 0
    tubeSectionDiameterCells =
        (tubeSectionDiameterCells == 0).select(1, tubeSectionDiameterCells);

    // Choose the best possible odd number from the Diameter array
    for (uint32_t diameterCounter = 1; diameterCounter <= numSections;
         ++diameterCounter) {
        // Verify if the cellsPerDiameter is odd or not
        if (tubeSectionDiameterCells(diameterCounter) % 2 == 0) {
            // Find the difference between the rounded and the actual
            // diameter value = Estimated diameter-Actual diameter
            const double diff =
                tubeSectionDiameterCells(diameterCounter) -
                (tubeSectionOriginalDiameter(diameterCounter) / ds);

            if (diff > 0) {
                tubeSectionDiameterCells(diameterCounter) -= 1;
            } else {
                tubeSectionDiameterCells(diameterCounter) += 1;
            }
        }
    }

    // Number of cells for total tube length
    const double   actualTubeLength = numSections * tube.segmentLength;
    const uint32_t totalTubeLengthInCells = std::round(actualTubeLength / ds);

    uint32_t domainX;
    uint32_t domainY;
    uint32_t domainZ;

    if (baffleSwitchFlag) {
        /// TODO: Implement the baffle condition;
        ///       For eccentric tube condition how should we include baffle?
    } else {
        domainX = totalTubeLengthInCells +
                  (1 + 1);  // +1 is excitation and +1 is dirichlet layer
        domainY = tubeSectionDiameterCells.maxCoeff() + 2;  // +2 is tube walls
        domainZ = tubeSectionDiameterCells.maxCoeff() + 2;  // +2 is tube walls
    }

    // Derive frame size
    uint32_t frameX;
    uint32_t frameY;
    uint32_t frameZ;

    frameX = domainX + 2;  // 2 = dead cells
    frameY = domainY + 2;  // 2 = dead cells
    frameZ = domainZ + 2;  // 2 = dead cells

    if (pmlSwitch) {
        frameX += 2 * pmlLayer;
        frameY += 2 * pmlLayer;
        frameZ += 2 * pmlLayer;
    }

    // Add all the frame params to the simGridParam
    SimulationGridParams simGridParams{frameX, frameY, frameZ};

    // Define PV_N array to store pressure and velocity components
    // PV_N(:,:,:, 1) = To store cell pressure
    // PV_N(:,:,:, 2) = To store Vx
    // PV_N(:,:,:, 3) = To store Vy
    // PV_N(:,:,:, 4) = To store Vz
    // PV_N(:,:,:, 5) = To store grid cell types
    auto PV_N = zeros<Tensor4<double>>(frameY, frameX, frameZ, 5);

    // Define cell types and store it in PV_N(:,:,:, 5)
    // Declare all the cells as air by default
    PV_N.chip(5, 3).setConstant(vt::cell_air);

    // Create the regular tube geometry

    // The axis of the cylindrical tube is along the x-axis, direction along
    // which acoustic wave properties. We will start building the tube model
    // from the starting position towards the end position.

    // Find the center of the tube along the yz-plane
    const uint32_t centerY = ceil(frameY / 2.0);
    const uint32_t centerZ = ceil(frameZ / 2.0);

    // Find the starting and ending point of the tube.
    // 1 is: for starting point, dead layer and excitation
    // 2 is: for cell_head in case of the circular baffle and to have an
    // empty cell between the cell_head and pml layers.
    StartInfo tubeStart;
    EndInfo   tubeEnd;

    if (true /* rad != 3 */) {
        tubeStart.startX = 1 + 1 + pmlLayer * pmlSwitch + 1;
        tubeEnd.endX = tubeStart.startX + totalTubeLengthInCells - 1;
    } else {
        /// TODO: check this later
    }

    tubeStart.startY = centerY;
    tubeEnd.endY = centerY;

    tubeStart.startZ = centerZ;
    tubeEnd.endZ = centerZ;

    // Store the cumulative length of each tube section
    auto tubeCummSectionLength = zeros<ArrayXd>(numSections);
    for (uint32_t sectionCount = 1; sectionCount <= numSections;
         ++sectionCount) {
        tubeCummSectionLength(sectionCount) = tube.segmentLength * sectionCount;
    }

    Array2X<uint32_t> currTubeSectionDiameterCells_SegmentCounter;
    if (junctionType == JunctionType::Centric) {
        currTubeSectionDiameterCells_SegmentCounter =
            vt::centricCircularContourGeneration(
                PV_N, numSections, totalTubeLengthInCells,
                tubeSectionDiameterCells, tubeCummSectionLength,
                boundaryInterpolation, tubeStart, ds);
    } else {
        currTubeSectionDiameterCells_SegmentCounter =
            vt::eccentricCircularContourGeneration(
                PV_N, numSections, totalTubeLengthInCells,
                tubeSectionDiameterCells, tubeCummSectionLength,
                boundaryInterpolation, tubeStart, pmlSwitch, pmlLayer, ds);
    }

    // Calculate tube midline length for eccentric geometry case
    // This can be done by checking the midpoint along y-axis
    double   vtMidlineLength = 0;
    uint32_t prevMidYPosition = 0;
    uint32_t minAirYInCells;
    uint32_t maxAirYInCells;

    for (uint32_t vtCrossSectionCounter = tubeStart.startX;
         vtCrossSectionCounter <= tubeEnd.endX; ++vtCrossSectionCounter) {
        gridPlaneProp = vt::findCellTypes(PV_N, vtCrossSectionCounter);

        // Find mid Y position
        minAirYInCells = UINT32_MAX;
        maxAirYInCells = 0;

        for (uint32_t y = 1; y <= rows(gridPlaneProp); ++y) {
            if (gridPlaneProp(y, tubeStart.startZ) ==
                GridCellTypeInplane::inVTContour) {
                if (y < minAirYInCells) {
                    minAirYInCells = y;
                }
                if (y > maxAirYInCells) {
                    maxAirYInCells = y;
                }
            }
        }

        const uint32_t currMidYPosition = (minAirYInCells + maxAirYInCells) / 2;

        if (vtCrossSectionCounter == tubeStart.startX) {
            prevMidYPosition = currMidYPosition;
            continue;
        }

        if (prevMidYPosition == currMidYPosition) {
            vtMidlineLength += ds;
            prevMidYPosition = currMidYPosition;
        } else {
            const double yPosDifferenceLen =
                fabs(currMidYPosition - prevMidYPosition) * ds;
            vtMidlineLength +=
                sqrt(ds * ds + yPosDifferenceLen * yPosDifferenceLen);
            prevMidYPosition = currMidYPosition;
        }
    }

    const double finalVtMidlineLength = vtMidlineLength + ds;

    // Percentage error in total tube length
    const double totalTubeLengthError =
        (finalVtMidlineLength - tube.totalLength) / tube.totalLength;
    fprintf(stderr, "Approximated tube length percentage error = %.4f \n",
            totalTubeLengthError * 100);

    // Define excitation cells
    // Check the grid cell types for the yz-plane that exists besides
    // (left-side) the tube starting point
    const uint32_t xExcitationPos = tubeStart.startX - 1;
    gridPlaneProp = vt::findCellTypes(PV_N, tubeStart.startX);

    // Find the grid size and traverse through yz-plane to assign
    // cell_excitation
    const auto gridSize = dimensions(PV_N);

    for (uint32_t yCount = 1; yCount <= gridSize[1]; ++yCount) {
        for (uint32_t zCount = 1; zCount <= gridSize[3]; ++zCount) {
            if (gridPlaneProp(yCount, zCount) ==
                GridCellTypeInplane::inVTContour) {
                PV_N(yCount, xExcitationPos, zCount, 5) = vt::cell_excitation;
            } else if (gridPlaneProp(yCount, zCount) ==
                       GridCellTypeInplane::onVTContour) {
                PV_N(yCount, xExcitationPos, zCount, 5) = vt::cell_wall;
            }
        }
    }

    // Define no_Pressure cells [Dirichlet Boundary Condition]
    if (rad == MouthTermination::DirichletBoundary) {
        tubeEnd.endX = tubeStart.startX + totalTubeLengthInCells - 1;
        gridPlaneProp = vt::findCellTypes(PV_N, tubeEnd.endX);

        // Traverse through yz-plane to assign cell_noPressure
        for (uint32_t yCount = 1; yCount <= gridSize[1]; ++yCount) {
            for (uint32_t zCount = 1; zCount <= gridSize[3]; ++zCount) {
                if (gridPlaneProp(yCount, zCount) ==
                    GridCellTypeInplane::inVTContour) {
                    PV_N(yCount, tubeEnd.endX + 1, zCount, 5) =
                        vt::cell_noPressure;
                } else if (gridPlaneProp(yCount, zCount) ==
                           GridCellTypeInplane::onVTContour) {
                    PV_N(yCount, tubeEnd.endX + 1, zCount, 5) =
                        vt::cell_noPressure;
                }
            }
        }
    }

    // Place two microphones:
    // 1. Near mouth-end [3mm inside of mouth] = listenerInfo
    // 2. Near glottal-end [Next to excitation cells] = sourceInfo
    // Deb: We need to careful for centric and eccentric tube.
    // First find the VT cross-section along the yz plane.
    // Then find the y coordinate for source and listener by checking
    // the min/max positions.

    // Determine the microphone/listener position near the mouth-end
    ListenerInfo listenerInfo;
    listenerInfo.listenerX = tubeEnd.endX - micXdistanceInCellsFromMouthEnd;
    listenerInfo.listenerZ = tubeEnd.endZ;
    gridPlaneProp = vt::findCellTypes(PV_N, listenerInfo.listenerX);

    // Find mid Y position
    minAirYInCells = UINT32_MAX;
    maxAirYInCells = 0;
    for (uint32_t y = 1; y <= rows(gridPlaneProp); ++y) {
        if (gridPlaneProp(y, tubeEnd.endZ) ==
            GridCellTypeInplane::inVTContour) {
            if (y < minAirYInCells) {
                minAirYInCells = y;
            }
            if (y > maxAirYInCells) {
                maxAirYInCells = y;
            }
        }
    }
    listenerInfo.listenerY = (minAirYInCells + maxAirYInCells) / 2;

    // Determine the source position near the glottal-end
    SourceInfo sourceInfo;
    sourceInfo.sourceX = tubeStart.startX + 1;
    sourceInfo.sourceZ = tubeStart.startZ;
    gridPlaneProp = vt::findCellTypes(PV_N, sourceInfo.sourceX);

    // Find mid Y position
    minAirYInCells = UINT32_MAX;
    maxAirYInCells = 0;
    for (uint32_t y = 1; y <= rows(gridPlaneProp); ++y) {
        if (gridPlaneProp(y, tubeStart.startZ) ==
            GridCellTypeInplane::inVTContour) {
            if (y < minAirYInCells) {
                minAirYInCells = y;
            }
            if (y > maxAirYInCells) {
                maxAirYInCells = y;
            }
        }
    }
    sourceInfo.sourceY = (minAirYInCells + maxAirYInCells) / 2;

    SimulationData data;
    data.gridParams = std::move(simGridParams);
    data.PV_N = std::move(PV_N);
    data.tubeStartArea = tubeStartArea;
    data.tubeStart = std::move(tubeStart);
    data.tubeEnd = std::move(tubeEnd);
    data.totalTubeLengthInCells = totalTubeLengthInCells;
    data.currTubeSectionDiameterCells_SegmentCounter =
        std::move(currTubeSectionDiameterCells_SegmentCounter);
    data.listenerInfo = std::move(listenerInfo);
    data.sourceInfo = std::move(sourceInfo);
    return data;
}
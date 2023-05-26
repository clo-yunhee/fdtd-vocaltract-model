#include "elliptical_tube_generation.hh"

#include "../constants.hh"
#include "area_function.hh"
#include "centric_elliptical_contour_generation.hh"
#include "eccentric_elliptical_contour_generation.hh"
#include "find_cell_types.hh"

vt::SimulationData vt::ellipticalTubeGeneration(
    SimulationType simulationType, JunctionType junctionType, Vowel vowel,
    bool boundaryInterpolation, bool baffleSwitchFlag, bool pmlSwitch,
    uint32_t pmlLayer, MouthTermination rad, double ds) {
    // Define tube geometry
    TubeProperties tube = vt::areaFunction(simulationType, junctionType, vowel);

    // Tube section area in m^2
    const auto tubeSectionArea_inm2 = tube.sectionArea_incm2 * 1e-2 * 1e-2;

    // Virtual mic/listener position near the mouth-end [inside the tube]
    const double   micXdistanceFromMouthEnd = 3e-3;  // 3 millimeters
    const uint32_t micXdistanceInCellsFromMouthEnd =
        (uint32_t)std::round(micXdistanceFromMouthEnd / ds);

    // Number of tube segments
    const uint32_t numSections = tubeSectionArea_inm2.elements();

    // Calculate tube cross-section semi-major and semi-minor axis length
    // The length of tube major and minor axis should be computed from the
    // original tube cross-sections/geometry [Arnela et al. 2016 - Influence
    // of vocal tract geometry simplifications]. However, in the case of
    // missing geometry images, the length of the semi-major and semi-minor
    // axis can be determined by defining a ratio by the user.

    bool isTubeGeometryDataExist = false;

    array majorAxisLen;
    array minorAxisLen;

    if (isTubeGeometryDataExist) {
        // TODO: We need vocal tract tube geometry / MRI images for this
    } else {
        // [Note: Keep semimajorAxisRatio > semiminorAxisRatio]
        const double semimajorAxisRatio = 3;
        const double semiminorAxisRatio = 1;

        // Determine length of the semi-major and semi-minor axis length
        // if semimajorAxis:semiminorAxis = a:b then their length can be
        // semimajorAxisLen = ax,  semiminorAxisLen = bx
        // ellipseArea = pi*ax*bx;
        // x = sqrt(ellipseArea/(semimajorAxisLen*semiminorAxisLen))

        auto lenX = sqrt(tubeSectionArea_inm2 /
                         (semimajorAxisRatio * semiminorAxisRatio * Pi));
        majorAxisLen = lenX * (semimajorAxisRatio * 2);
        minorAxisLen = lenX * (semiminorAxisRatio * 2);
    }

    array tubeMajorAxisLenInCells = round(majorAxisLen / ds).as(u32);
    array tubeMinorAxisLenInCells = round(minorAxisLen / ds).as(u32);

    // Cross-sectional area at the tube start end
    // I didn't use "tubeSectionArea_inm2" directly because I wanted to use
    // the derived area by using the derived diameter.
    const double tubeStartArea =
        (Pi * (majorAxisLen(0) / 2) * (minorAxisLen(0) / 2)).scalar<double>();

    // Change the majorAxisLenInCells and minorAxisLenInCells to 1
    // if it contains 0
    tubeMajorAxisLenInCells(tubeMajorAxisLenInCells == 0) = 1;
    tubeMinorAxisLenInCells(tubeMinorAxisLenInCells == 0) = 1;

    // Choose the best possible odd number from the MinorAxisLenInCell
    // and tubeMajorAxisLenInCells array. In this model, the major axis goes
    // along the xz-plane.
    for (uint32_t counter = 0; counter < numSections; ++counter) {
        // Verify if the tubeMajorAxisLenInCells is odd or not
        if (allTrue<bool>(tubeMajorAxisLenInCells(counter) % 2 == 0)) {
            // Find the difference between the rounded and the actual
            // diameter value = Estimated diameter-Actual diameter
            const double diff = (tubeMajorAxisLenInCells(counter) -
                                 (majorAxisLen(counter) / ds))
                                    .scalar<double>();

            if (diff > 0) {
                tubeMajorAxisLenInCells(counter) =
                    tubeMajorAxisLenInCells(counter) - 1;
            } else {
                tubeMajorAxisLenInCells(counter) =
                    tubeMajorAxisLenInCells(counter) + 1;
            }
        }
    }

    // Store the ellipse semi-major axis and semi-minor axis diameter
    array ellipseAxisLenInfo = constant(0, 2, numSections, u32);
    ellipseAxisLenInfo(0, span) = tubeMajorAxisLenInCells;
    ellipseAxisLenInfo(1, span) = tubeMinorAxisLenInCells;

    // Find the total tube length and calculate the percentage error
    // in the approximated tube length

    // Number of cells for total tube length
    const double   actualTubeLength = numSections * tube.segmentLength;
    const uint32_t totalTubeLengthInCells =
        (uint32_t)std::round(actualTubeLength / ds);

    uint32_t domainX;
    uint32_t domainY;
    uint32_t domainZ;

    if (baffleSwitchFlag) {
        /// TODO: Implement the baffle condition;
        ///       For eccentric tube condition how should we include baffle?
    } else {
        // The major axis of the elliptical cross-section goes along the
        // z-axis and the minor axis of the elliptical cross-section goes
        // along the y-axis.
        domainX = totalTubeLengthInCells +
                  (1 + 1);  // +1 is excitation and +1 is dirichlet layer
        domainY =
            max<uint32_t>(tubeMinorAxisLenInCells) + 2;  // +2 is tube walls
        domainZ =
            max<uint32_t>(tubeMajorAxisLenInCells) + 2;  // +2 is tube walls
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
    array PV_N = constant(0, frameY, frameX, frameZ, 5);

    // Define cell types and store it in PV_N(:,:,:, 5)
    // Declare all the cells as air by default
    PV_N(span, span, span, 4) = vt::cell_air;

    // Create the regular tube geometry

    // The axis of the cylindrical tube is along the x-axis, direction along
    // which acoustic wave properties. We will start building the tube model
    // from the starting position towards the end position.

    // Find the center of the tube along the yz-plane
    const uint32_t centerY = (uint32_t)ceil((frameY - 1) / 2.0);
    const uint32_t centerZ = (uint32_t)ceil((frameZ - 1) / 2.0);

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
    array tubeCummSectionLength = constant(0, 1, numSections);
    for (uint32_t sectionCount = 0; sectionCount < numSections;
         ++sectionCount) {
        tubeCummSectionLength(sectionCount) =
            tube.segmentLength * (sectionCount + 1);
    }

    array currTubeSectionDiameterCells_SegmentCounter;
    if (junctionType == JunctionType::Centric) {
        currTubeSectionDiameterCells_SegmentCounter =
            vt::centricEllipticalContourGeneration(
                PV_N, numSections, totalTubeLengthInCells,
                tubeCummSectionLength, ellipseAxisLenInfo,
                boundaryInterpolation, tubeStart, ds);
    } else {
        currTubeSectionDiameterCells_SegmentCounter =
            vt::eccentricEllipticalContourGeneration(
                PV_N, numSections, totalTubeLengthInCells,
                tubeCummSectionLength, ellipseAxisLenInfo,
                boundaryInterpolation, tubeStart, pmlSwitch, pmlLayer, ds);
    }

    // Calculate tube midline length for eccentric geometry case
    // This can be done by checking the midpoint along y-axis
    double   vtMidlineLength = 0;
    uint32_t prevMidYPosition = 0;

    for (uint32_t vtCrossSectionCounter = tubeStart.startX;
         vtCrossSectionCounter <= tubeEnd.endX; ++vtCrossSectionCounter) {
        const array gridPlaneProp =
            vt::findCellTypes(PV_N, vtCrossSectionCounter);
        const array    airCells = where(gridPlaneProp(span, tubeStart.startZ) ==
                                        (float)GridCellTypeInplane::inVTContour);
        const uint32_t currMidYPosition =
            ((min(airCells) + max(airCells)) / 2).scalar<uint32_t>();

        if (vtCrossSectionCounter == tubeStart.startX) {
            prevMidYPosition = currMidYPosition;
            continue;
        }

        if (prevMidYPosition == currMidYPosition) {
            vtMidlineLength += ds;
            prevMidYPosition = currMidYPosition;
        } else {
            const double yPosDifferenceLen =
                std::abs((double)currMidYPosition - (double)prevMidYPosition) *
                ds;
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
    const array    gridPlaneProp = vt::findCellTypes(PV_N, tubeStart.startX);

    // Find the grid size and traverse through yz-plane to assign
    // cell_excitation
    const dim4 gridSize = PV_N.dims();

    for (uint32_t yCount = 0; yCount < gridSize[0]; ++yCount) {
        for (uint32_t zCount = 0; zCount < gridSize[2]; ++zCount) {
            if (allTrue<bool>(gridPlaneProp(yCount, zCount) ==
                              (float)GridCellTypeInplane::inVTContour)) {
                PV_N(yCount, xExcitationPos, zCount, 4) = vt::cell_excitation;
            } else if (allTrue<bool>(gridPlaneProp(yCount, zCount) ==
                                     (float)GridCellTypeInplane::onVTContour)) {
                PV_N(yCount, xExcitationPos, zCount, 4) = vt::cell_wall;
            }
        }
    }

    // Define no_Pressure cells [Dirichlet Boundary Condition]
    if (rad == MouthTermination::DirichletBoundary) {
        tubeEnd.endX = tubeStart.startX + totalTubeLengthInCells - 1;
        const array gridPlaneProp = vt::findCellTypes(PV_N, tubeEnd.endX);

        // Traverse through yz-plane to assign cell_noPressure
        for (uint32_t yCount = 0; yCount < gridSize[0]; ++yCount) {
            for (uint32_t zCount = 0; zCount < gridSize[2]; ++zCount) {
                if (allTrue<bool>(gridPlaneProp(yCount, zCount) ==
                                  (float)GridCellTypeInplane::inVTContour)) {
                    PV_N(yCount, tubeEnd.endX + 1, zCount, 4) =
                        vt::cell_noPressure;
                } else if (allTrue<bool>(
                               gridPlaneProp(yCount, zCount) ==
                               (float)GridCellTypeInplane::onVTContour)) {
                    PV_N(yCount, tubeEnd.endX + 1, zCount, 4) =
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
    {
        listenerInfo.listenerX = tubeEnd.endX - micXdistanceInCellsFromMouthEnd;
        listenerInfo.listenerZ = tubeEnd.endZ;
        const array gridPlaneProp =
            vt::findCellTypes(PV_N, listenerInfo.listenerX);
        const array airCells =
            where(gridPlaneProp(span, listenerInfo.listenerZ) ==
                  (float)GridCellTypeInplane::inVTContour);
        const uint32_t yPosition =
            ((min(airCells) + max(airCells)) / 2).scalar<uint32_t>();
        listenerInfo.listenerY = yPosition;
    }

    // Determine the source position near the glottal-end
    SourceInfo sourceInfo;
    {
        sourceInfo.sourceX = tubeStart.startX + 1;
        sourceInfo.sourceZ = tubeStart.startZ;
        const array gridPlaneProp = vt::findCellTypes(PV_N, sourceInfo.sourceX);
        const array airCells = where(gridPlaneProp(span, sourceInfo.sourceZ) ==
                                     (float)GridCellTypeInplane::inVTContour);
        const uint32_t yPosition =
            ((min(airCells) + max(airCells)) / 2).scalar<uint32_t>();
        sourceInfo.sourceY = yPosition;
    }

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
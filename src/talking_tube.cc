#include <af/util.h>

#include "audio/generate_audio.hh"
#include "audio/save_audio_file.hh"
#include "constants.hh"
#include "image/make_cell_image.hh"
#include "image/make_wave_image.hh"
#include "image/save_to_png.hh"
#include "model/LF.h"
#include "routines/circular_tube_generation.hh"
#include "routines/elliptical_tube_generation.hh"
#include "routines/impulse_signal.hh"
#include "types.hh"

using namespace vt;

void talkingTube() {
    const int              srateBase = 44100;
    const uint32_t         srateMul = 12;  // Sample rate multiplier
    const double           dur = 50e-3;    // simulation duration
    const SimulationType   simulationType = SimulationType::VowelSound;
    const Vowel            vowel = Vowel::I;
    const CrossSectionType crossSectionType = CrossSectionType::Circular;
    const JunctionType     junctionType = JunctionType::Centric;
    const SourceType       sourceModelType = SourceType::GFM_LF;
    const double           LF_Rd = 2.4;
    const bool             pmlSwitch = true;
    const bool             baffleSwitchFlag = false;
    const MouthTermination rad = MouthTermination::DirichletBoundary;
    const bool             boundaryInterpolation = false;
    const bool             plotCells = true;
    const bool             plotPressure = true;

    // DEFINE CONSTANTS

    const double maxSigmadt = 0.5;  // Attenuation coefficient at the PML layer
    const double srate = srateBase * srateMul;  // Sample frequency
    const uint32_t pmlLayer = 6;                // Number of PML layers

    // DASHBOARD

    const double dt = 1 / srate;  // Temporal resolution / sample time period
    const double dx = dt * c * sqrt_3;  // Spatial resolution along x-direction
    const double dy = dt * c * sqrt_3;  // Spatial resolution along y-direction
    const double dz = dt * c * sqrt_3;  // Spatial resolution along z-direction
    const double AudioTime = dur;
    const double kappa = rho * c * c;  // Bulk modulus
    const double ds = dx;              // Spatial resoltuion(ds) = dx = dy
    const double rhoC = rho * c;

    // SIMULATION TIME
    std::vector<double> times;
    uint32_t            timeIndex = 0;
    double              curTime = 0;
    while (curTime < AudioTime - dt) {
        times.push_back(curTime);
        timeIndex++;
        curTime = timeIndex * dt;
    }
    const uint32_t STEPS = times.size();
    array          t(1, STEPS, times.data());

    // DEFINE BETA AND SIGMAPRIME PARAMETERS

    double vis_Min, vis_Max, vis_Boundary;

    if (sourceModelType == SourceType::Impulse) {
        vis_Min = -10000;
        vis_Max = 40000;
        vis_Boundary = 40000;
    } else if (sourceModelType == SourceType::Sine) {
        vis_Min = -4000;
        vis_Max = 4000;
        vis_Boundary = 4000;
    } else if (sourceModelType == SourceType::GFM_LF) {
        vis_Min = -2000;
        vis_Max = 25000;
        vis_Boundary = 25000;
    }

    array sigmadt = constant(0, pmlLayer);

    /*
     To store beta(tube wall) and sigmaPrimedt(PML Layers) [beta, sigmaPrimedt].
     beta_air = 1, beta_PML = 1 and beta_wall = 0
     sigmaPrimedt = sigmaPrime*dt
     sigmaPrime = 1 - beta + sigma
     e.g -
     sigma=0 for all the non-PML layers. Hence, sigmaPrime = 1 - beta
     inside the domain. Therefore,
     WALL -> beta = 0, sigma_prima*dt = (1-0)*dt = 1*dt = dt
     AIR  -> beta = 1, sigma_prima*dt = (1-1)*dt = 0*dt = 0
     HEAD CELLS -> beta = 0, sigma_prima*dt = (1-0)*dt = 1*dt = dt
     NOTE - We are considering excitation cells/head cells as special wall cells
    */
    array typeValues = constant(0, 2, vt::cell_numTypes);
    typeValues(span, vt::cell_wall) = array(dim4(2, 1), {0., dt});  // VT walls
    typeValues(span, vt::cell_air) = array(dim4(2, 1), {1, 0});     // air
    typeValues(span, vt::cell_noPressure) = array(dim4(2, 1), {1, 0});  // air
    typeValues(span, vt::cell_excitation) =
        array(dim4(2, 1), {0., dt});  // excitation
    typeValues(span, vt::cell_head) =
        array(dim4(2, 1), {0., dt});  // head cells

    // Define beta and sigmaPrimedt for PML layers
    // For PML layers beta = 1,
    // sigmaPrimedt = [1-beta+sigma]*dt = [1-1+sigma]*dt = sigmadt
    for (uint32_t pmlCounter = 0; pmlCounter <= pmlLayer - 1; ++pmlCounter) {
        sigmadt(pmlCounter) =
            ((double)pmlCounter / (pmlLayer - 1)) * maxSigmadt;
        typeValues(span, vt::cell_pml0 + pmlCounter) = array(
            dim4(2, 1), {1., (double)sigmadt(pmlCounter).scalar<float>()});
    }
    typeValues(span, vt::cell_dead) = array(dim4(2, 1), {0, 1000000});

    // DEFINE SIMULATION TYPE

    // [Note]: I have created two separate functions for the circular and
    // elliptical tube generation. However, that's not required since circle
    // is a special case of ellipse where the ratio between the
    // semiMajorAxis and semiMinorAxis is equal.
    // [semiMajorAxis:semiMinorAxis = 1:1]

    SimulationData simData;

    switch (simulationType) {
        case SimulationType::OpenSpace:
            break;
        //
        case SimulationType::RegularTube:
        case SimulationType::VowelSound:
            switch (crossSectionType) {
                case CrossSectionType::Circular:
                    simData = vt::circularTubeGeneration(
                        simulationType, junctionType, vowel,
                        boundaryInterpolation, baffleSwitchFlag, pmlSwitch,
                        pmlLayer, rad, ds);
                    break;
                case CrossSectionType::Elliptical:
                    simData = vt::ellipticalTubeGeneration(
                        simulationType, junctionType, vowel,
                        boundaryInterpolation, baffleSwitchFlag, pmlSwitch,
                        pmlLayer, rad, ds);
                    break;
                case CrossSectionType::Square:
                    /*simData = vt::squareTubeGeneration(
                        simulationType, junctionType, vowel,
                        boundaryInterpolation, baffleSwitchFlag, pmlSwitch,
                        pmlLayer, rad, ds);*/
                    throw std::invalid_argument(
                        "square cross section not implemented yet");
                    break;
                default:
                    throw std::invalid_argument("unknown cross section type");
            }
            break;
        default:
            throw std::invalid_argument("unknown simulation type");
    }

    // VARIABLES TO SAVE PRESSURE AT MOUTH_END AND GLOTTAL_END

    array Pr_mouthEnd = constant(0, 1, STEPS);
    array Pr_glottalEnd = constant(0, 1, STEPS);

    // DETERMINE WALL IMPEDANCE

    const double mu3D = 0.005;  // Boundary Admittance Coefficient
    const double alpha3D =
        1 / (0.5 + 0.25 * (mu3D + (1 / mu3D)));  // Sound Absorption Coefficient
    const double z_inv =
        1 / (rhoC * ((1 + sqrt(1 - alpha3D)) / (1 - sqrt(1 - alpha3D))));

    // SOURCE MODEL

    array excitationV = constant(0, 1, STEPS);

    switch (sourceModelType) {
        case SourceType::Sine: {
            // Sine wave source model
            const double excitationF = 440;  // 440 Hz
            const double srcAmplitude = 25;
            for (uint32_t i = 0; i < STEPS; ++i) {
                excitationV(i) =
                    srcAmplitude * sin(2 * Pi * excitationF * dt * i);
            }
        } break;
        case SourceType::Gaussian: {
            // Gaussian source model
            const double f0 = 10e3;  // 10 kHz
            const double bellPeakPos = 0.646 / f0;
            const double bellWidth = 0.29 * bellPeakPos;
            excitationV = exp(-pow((t - bellPeakPos) / bellWidth, 2.));
        } break;
        case SourceType::Impulse: {
            // Impulse response function
            excitationV = src::ImpulseSignal(srate, 10000, 2, srateBase / 2);
        } break;
        case SourceType::GaussianNoise: {
            // Gaussian white noise
            // excitationV = randn(STEPS, 1);
        } break;
        case SourceType::VFModel: {
            // Vocal fold model - Two Mass Model
            // Set the vocal fold parameters
            /* [airParam, vf_structuralParam, vf_flowParam, vf_matParam] =
               vf_SetVocalFoldParams(); */
        } break;
        case SourceType::GFM_LF: {
            // Glottal flow model - LF model

            const double excitationF = 230;
            const double srcAmplitude = 4000;

            models::LF model;
            model.setRd(LF_Rd);

            for (uint32_t i = 0; i < STEPS; ++i) {
                excitationV(i) = srcAmplitude * model.evaluateAntiderivative(
                                                    excitationF * dt * i);
            }
        } break;
        default:
            throw std::invalid_argument("unknown source model type");
    }

    // Define source propagation direction
    // srcDirection index: 1 = Left  = -1
    //                     2 = Down  = -1
    //                     3 = Back  = -1
    //                     4 = Right =  1
    //                     5 = Up    =  1
    //                     6 = Front =  1
    array srcDirection =
        array(dim4(6), {0, 0, 0, 1, 0, 0});  // For forward direction

    // SOURCE AND VIRTUAL MIC

    const uint32_t frameX = simData.gridParams.frameX;
    const uint32_t frameY = simData.gridParams.frameY;
    const uint32_t frameZ = simData.gridParams.frameZ;

    // Find the mid point of the domain for open-space simulation
    const uint32_t midX = frameX / 2;
    const uint32_t midY = frameY / 2;
    const uint32_t midZ = frameZ / 2;

    // DEFINE PML LAYERS AND DEAD CELLS INSIDE THE 3D GRID

    // Define cell_dead to the outer-most layer of the frame

    // ---Adding dead cells to the top and bottom surfaces---
    simData.PV_N(0, seq(0, frameX - 1), seq(0, frameZ - 1), 4) = vt::cell_dead;
    simData.PV_N(frameY - 1, seq(0, frameX - 1), seq(0, frameZ - 1), 4) =
        vt::cell_dead;

    // ---Adding dead cells to the back and front surfaces---
    simData.PV_N(seq(0, frameY - 1), seq(0, frameX - 1), 0, 4) = vt::cell_dead;
    simData.PV_N(seq(0, frameY - 1), seq(0, frameX - 1), frameZ - 1, 4) =
        vt::cell_dead;

    // ---Adding dead cells to the left and right surfaces---
    simData.PV_N(seq(0, frameY - 1), 0, seq(0, frameZ - 1), 4) = vt::cell_dead;
    simData.PV_N(seq(0, frameY - 1), frameX - 1, seq(0, frameZ - 1), 4) =
        vt::cell_dead;

    if (pmlSwitch) {
        // Define PML layers starting from the outer-most layer
        uint32_t pmlType = (uint32_t)vt::cell_pml5;

        const uint32_t xShift = 1;
        uint32_t       xStart = 1;
        uint32_t       xEnd = frameX - 2;

        const uint32_t yShift = 1;
        uint32_t       yStart = 1;
        uint32_t       yEnd = frameY - 2;

        const uint32_t zShift = 1;
        uint32_t       zStart = 1;
        uint32_t       zEnd = frameZ - 2;

        for (uint32_t pmlCount = 1; pmlCount <= pmlLayer; ++pmlCount) {
            // ---Adding pml layers to the top and bottom surfaces---
            simData.PV_N(yShift + pmlCount, seq(xStart, xEnd),
                         seq(zStart, zEnd), 4) = pmlType;
            simData.PV_N(frameY - 1 - pmlCount, seq(xStart, xEnd),
                         seq(zStart, zEnd), 4) = pmlType;

            // ---Adding pml layers to the back and front surfaces---
            simData.PV_N(seq(yStart, yEnd), seq(xStart, xEnd),
                         zShift + pmlCount, 4) = pmlType;
            simData.PV_N(seq(yStart, yEnd), seq(xStart, xEnd),
                         frameZ - 1 - pmlCount, 4) = pmlType;

            // ---Adding pml layers to the left and right surfaces---
            simData.PV_N(seq(yStart, yEnd), xShift + pmlCount,
                         seq(zStart, zEnd), 4) = pmlType;
            simData.PV_N(seq(yStart, yEnd), frameX - 1 - pmlCount,
                         seq(zStart, zEnd), 4) = pmlType;

            // Update the pmlLayer type.
            // Update the start and end position
            pmlType = pmlType - 1;
            xStart = xStart + 1;
            xEnd = xEnd - 1;

            yStart = yStart + 1;
            yEnd = yEnd - 1;

            zStart = zStart + 1;
            zEnd = zEnd - 1;
        }
    }

    // MODEL VISUALISATION

    if (plotCells) {
        // Visualize the PML layers and excitation cell by taking a slice of the
        // domain
        array buildFrame = simData.PV_N(span, span, span, 4);
        // Uncomment to visualise source position for the open space simulation
        // buildFrame(midY, midX, midZ) = -1;
        buildFrame(simData.listenerInfo.listenerY,
                   simData.listenerInfo.listenerX,
                   simData.listenerInfo.listenerZ) = GridCellType::cell_head;

        // Slice along XZ-plane
        array frame_xz =
            buildFrame((uint32_t)std::ceil((frameY - 1) / 2.0), span, span);

        // Slice along YZ-plane
        array frame_yz =
            buildFrame(span, (uint32_t)std::ceil((frameX - 1) / 2.0), span);

        // Slice along XY-plane
        array frame_xy =
            buildFrame(span, span, (uint32_t)std::ceil((frameZ - 1) / 2.0));

        // Save images
        const auto frame_xz_image =
            image::makeCellImage(moddims(frame_xz, frameX, frameZ));
        const auto frame_yz_image =
            image::makeCellImage(moddims(frame_yz, frameY, frameZ));
        const auto frame_xy_image =
            image::makeCellImage(moddims(frame_xy, frameY, frameX));
        image::saveToPng("images/frame_xz.png", frame_xz_image);
        image::saveToPng("images/frame_yz.png", frame_yz_image);
        image::saveToPng("images/frame_xy.png", frame_xy_image);
    }

    // LOOK UP TABLE FOR BOUNDARY CONDITIONS

    // Look up table store cell types in the following order:
    // [right_current, up_current, front_current]
    const GridCellType lookUp_table[][3] = {
        {cell_air, cell_air, cell_air},   {cell_air, cell_air, cell_wall},
        {cell_air, cell_wall, cell_air},  {cell_air, cell_wall, cell_wall},
        {cell_wall, cell_air, cell_air},  {cell_wall, cell_air, cell_wall},
        {cell_wall, cell_wall, cell_air}, {cell_wall, cell_wall, cell_wall},
    };
    const uint32_t lookUp_table_count = 8;

    const double air_normalV_component[] = {0, 1,      1,      0.7071,
                                            1, 0.7071, 0.7071, 0.5774};

    const double wall_normalV_component[] = {0.5774, 1,      1,      0.7071,
                                             1,      0.7071, 0.7071, 0};

    // CALCULATE MINBETA AND MAXSIGMAPRIMEDT
    // -> Store the cellTypes, typeIndex, beta, and sigmaPrimedt separately
    // -> Compute minbeta and maxsigmaprimedt

    const uint32_t height =
        frameY - 2;  // => Number of cells along Y - cell_Dead
    const uint32_t width =
        frameX - 2;  // => Number of cells along X - cell_Dead
    const uint32_t depth =
        frameZ - 2;  // => Number of cells along Z - cell_Dead

    // Create arrays to store cellTypes, beta, and sigmaPrime_dt values
    // We create 4 columns here as each grid cell has three neighbours
    // Note => We don't need to define these arrays [Just for my understanding]
    array cellTypes = constant(0, 1, 4);
    array beta = constant(0, 1, 4);
    array sigmaPrime_dt = constant(0, 1, 4);

    // To store is_excitation and excitation_weight
    const uint32_t numGridCellsForComputation =
        (frameX - 2) * (frameY - 2) * (frameZ - 2);

    array is_excitation = constant(0, 4, numGridCellsForComputation);
    array excitation_weight = constant(0, 3, numGridCellsForComputation);
    array xor_term = constant(0, 6, numGridCellsForComputation);
    array N_out = constant(0, 1, numGridCellsForComputation);
    array N_in = constant(0, 1, numGridCellsForComputation);

    // Create arrays to store minBeta and maxSigmaPrime_dt
    // Note => We consider minBeta, because beta parameter is defined at the
    // centre of the cell. And the velocity components are defined either on
    // edges [2D grid cell]/ side surfaces [3D grid cell]. Therefore, we
    // tie-break between these grid cells considering minBeta value.

    array minVxBeta = constant(0, height, width, depth);
    array minVyBeta = constant(0, height, width, depth);
    array minVzBeta = constant(0, height, width, depth);

    array maxVxSigmaPrimedt = constant(0, height, width, depth);
    array maxVySigmaPrimedt = constant(0, height, width, depth);
    array maxVzSigmaPrimedt = constant(0, height, width, depth);

    // Store the excitation_weight and are_we_not_excitations in matrix format
    // to compute Vx, Vy, Vz
    array excitation_Vx_weight = constant(0, height, width, depth);
    array excitation_Vy_weight = constant(0, height, width, depth);
    array excitation_Vz_weight = constant(0, height, width, depth);

    array are_we_not_excitations_Vx = constant(0, height, width, depth);
    array are_we_not_excitations_Vy = constant(0, height, width, depth);
    array are_we_not_excitations_Vz = constant(0, height, width, depth);

    // Store the xor_val in matrix format
    array xor_val1 = constant(0, height, width, depth);
    array xor_val2 = constant(0, height, width, depth);
    array xor_val3 = constant(0, height, width, depth);
    array xor_val4 = constant(0, height, width, depth);
    array xor_val5 = constant(0, height, width, depth);
    array xor_val6 = constant(0, height, width, depth);

    // Store the N_out and N_in in matrix format
    array N_out_mat = constant(0, height, width, depth);
    array N_in_mat = constant(0, height, width, depth);

    // Note => We are reshaping sigmaPrimedt in matrix format as we need this
    // while calculating pressure [Check the denominator of the discretized
    // pressure equation]
    array sigmaPrimedt = constant(0, height, width, depth);

    // Note => Here for each grid cell, we store cell type, beta and
    // siggmaPrine_dt values in arrays for its neighbor cells and the current
    // cell. In case of For example, in case of -
    // 2D/2.5D => [current, right_current, top_current]
    // 3D => [current, right_current, up_current, front_current]

    // Set a counter to save data for is_excitation, excitation_weight, xor,
    // N_out, N_in
    uint32_t counter = 0;

    for (uint32_t height_idx = 1; height_idx <= frameY - 2; ++height_idx) {
        for (uint32_t width_idx = 1; width_idx <= frameX - 2; ++width_idx) {
            for (uint32_t depth_idx = 1; depth_idx <= frameZ - 2; ++depth_idx) {
                // Find the cellTypes = [current, right_current, up_current,
                // front_current]
                cellTypes(0) =
                    simData.PV_N(height_idx, width_idx, depth_idx, 4);
                cellTypes(1) =
                    simData.PV_N(height_idx, width_idx + 1, depth_idx, 4);
                cellTypes(2) =
                    simData.PV_N(height_idx - 1, width_idx, depth_idx, 4);
                cellTypes(3) =
                    simData.PV_N(height_idx, width_idx, depth_idx + 1, 4);

                // Store beta values in beta array
                beta = typeValues(0, cellTypes);

                // Calculate minBeta
                const array min_beta_Vx =
                    min(beta(0), beta(1));  // minBeta(current, right_current)
                const array min_beta_Vy =
                    min(beta(0), beta(2));  // minBeta(current, top_current)
                const array min_beta_Vz =
                    min(beta(0), beta(3));  // minBeta(current, front_current)

                minVxBeta(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    min_beta_Vx;
                minVyBeta(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    min_beta_Vy;
                minVzBeta(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    min_beta_Vz;

                // Store sigmaPrime*dt values in sigmaPrime_dt
                sigmaPrime_dt = typeValues(1, cellTypes);

                // Store sigmaPrime_dt of the current cell only
                sigmaPrimedt(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    sigmaPrime_dt(0);

                // Calculate maxSigmaPrime_dt
                const array max_sigmaPrimedt_Vx =
                    max(sigmaPrime_dt(0), sigmaPrime_dt(1));
                const array max_sigmaPrimedt_Vy =
                    max(sigmaPrime_dt(0), sigmaPrime_dt(2));
                const array max_sigmaPrimedt_Vz =
                    max(sigmaPrime_dt(0), sigmaPrime_dt(3));

                maxVxSigmaPrimedt(height_idx - 1, width_idx - 1,
                                  depth_idx - 1) = max_sigmaPrimedt_Vx;
                maxVySigmaPrimedt(height_idx - 1, width_idx - 1,
                                  depth_idx - 1) = max_sigmaPrimedt_Vy;
                maxVzSigmaPrimedt(height_idx - 1, width_idx - 1,
                                  depth_idx - 1) = max_sigmaPrimedt_Vz;

                // Check whether the current cell is an excitation cell or not
                is_excitation(span, counter) =
                    (cellTypes == (uint32_t)vt::cell_excitation);

                // Verify if both the current cell and neighbouring cells are
                // not excitation cells
                are_we_not_excitations_Vx(height_idx - 1, width_idx - 1,
                                          depth_idx - 1) =
                    (1 - is_excitation(0, counter)) *
                    (1 - is_excitation(1, counter));

                are_we_not_excitations_Vy(height_idx - 1, width_idx - 1,
                                          depth_idx - 1) =
                    (1 - is_excitation(0, counter)) *
                    (1 - is_excitation(2, counter));

                are_we_not_excitations_Vz(height_idx - 1, width_idx - 1,
                                          depth_idx - 1) =
                    (1 - is_excitation(0, counter)) *
                    (1 - is_excitation(3, counter));

                // Determine excitation weight for the current cell
                // Find excitation velocity going out of the cell and coming
                // back into the cells depending upon the source direction.
                // Add them all to find the net velocity [velocity components]
                // for the current cell.

                const array excitationWeightForward =
                    tile(is_excitation(0, counter), 3) *
                    srcDirection(seq(3, 5));

                const array excitationWeightBackward =
                    is_excitation(seq(1, 3), counter) * srcDirection(seq(0, 2));

                excitation_weight(span, counter) =
                    excitationWeightForward + excitationWeightBackward;

                // Store the excitation_weight in matrix format
                excitation_Vx_weight(height_idx - 1, width_idx - 1,
                                     depth_idx - 1) =
                    excitation_weight(0, counter);

                excitation_Vy_weight(height_idx - 1, width_idx - 1,
                                     depth_idx - 1) =
                    excitation_weight(1, counter);

                excitation_Vz_weight(height_idx - 1, width_idx - 1,
                                     depth_idx - 1) =
                    excitation_weight(2, counter);

                // Check the adjacent cells to determine the velocity direction
                xor_val1(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(1) * (1 - beta(0));
                xor_val2(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(0) * (1 - beta(1));

                xor_val3(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(2) * (1 - beta(0));
                xor_val4(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(0) * (1 - beta(2));

                xor_val5(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(3) * (1 - beta(0));
                xor_val6(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(0) * (1 - beta(3));

                // Find the boundary condition from the lookup table and
                // select the corresponding unit vector
                uint32_t lookup_idx = (uint32_t)(-1);

                for (uint32_t i = 0; i < lookUp_table_count; ++i) {
                    if (allTrue<bool>(lookUp_table[i][0] == beta(1) &&
                                      lookUp_table[i][1] == beta(2) &&
                                      lookUp_table[i][2] == beta(3))) {
                        lookup_idx = i;
                        break;
                    }
                }

                if (lookup_idx == (uint32_t)(-1)) {
                    throw std::logic_error("invalid state");
                }

                N_out(0, counter) = air_normalV_component[lookup_idx] * beta(0);
                N_in(0, counter) =
                    wall_normalV_component[lookup_idx] * (1 - beta(0));

                N_out_mat(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    N_out(0, counter);
                N_in_mat(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    N_in(0, counter);

                // Increment the counter
                counter++;
            }
        }
    }

    const array betaVxSqr = pow(minVxBeta, 2);
    const array betaVySqr = pow(minVyBeta, 2);
    const array betaVzSqr = pow(minVzBeta, 2);

    const array betaVxSqr_dt_invRho = (betaVxSqr * dt) / rho;
    const array betaVySqr_dt_invRho = (betaVySqr * dt) / rho;
    const array betaVzSqr_dt_invRho = (betaVzSqr * dt) / rho;

    const double rho_sqrC_dt_invds = (kappa * dt) / ds;
    const double rho_sqrC_dt = kappa * dt;

    // INITIALISE ACOUSTIC PARAMETERS AND COEFFICIENTS
    array PV_NPlus1 = constant(0, frameY, frameX, frameZ, 5);

    array Ug_array = constant(0, 1, STEPS);

    array CxVx = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array CyVy = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array CzVz = constant(0, frameY - 2, frameX - 2, frameZ - 2);

    array CxP = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array CyP = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array CzP = constant(0, frameY - 2, frameX - 2, frameZ - 2);

    array Pr_next = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array Vx_next = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array Vy_next = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array Vz_next = constant(0, frameY - 2, frameX - 2, frameZ - 2);

    array vb_alphaX = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array vb_alphaY = constant(0, frameY - 2, frameX - 2, frameZ - 2);
    array vb_alphaZ = constant(0, frameY - 2, frameX - 2, frameZ - 2);

    // Introduce Dirichlet boundary condition
    uint32_t xNoPressure;
    array    isNoPressureCell;

    if (rad == MouthTermination::DirichletBoundary &&
        simulationType != SimulationType::OpenSpace) {
        xNoPressure = simData.tubeStart.startX + simData.totalTubeLengthInCells;
        // gridNoPressurePlane = PV_N(:, xNoPressure, :, 5)
        const array gridNoPressurePlane =
            simData.PV_N(span, xNoPressure, span, 4);
        isNoPressureCell =
            (gridNoPressurePlane == (float)GridCellType::cell_noPressure);
    }

    // MOVE TO GPU

    // --

    for (uint32_t T = 0; T < STEPS; ++T) {
        // STEP1: Calculate [del.V] = [dVx/dx + dVy/dy + dVz/dz]
        // CxVx = dVx, where Vx = velocity along x direction = Vx_curr - Vx_left
        // CyVy = dVy, where Vy = velocity along y direction = Vy_curr - Vy_down
        // CzVz = dVz, where Vz = velocity along z direction = Vz_curr - Vz_back

        CxVx(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 1) -
            simData.PV_N(seq(1, frameY - 2), seq(0, frameX - 3),
                         seq(1, frameZ - 2), 1);

        CyVy(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 2) -
            simData.PV_N(seq(2, frameY - 1), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 2);

        CzVz(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 3) -
            simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(0, frameZ - 3), 3);

        // STEP2: Calculate Pr_next

        Pr_next(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            (simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                          seq(1, frameZ - 2), 0) -
             (rho_sqrC_dt_invds * (CxVx + CyVy + CzVz))) /
            (1 + sigmaPrimedt);

        // Make pressure values of no_pressure cells as zeros
        if (rad == MouthTermination::DirichletBoundary &&
            simulationType != SimulationType::OpenSpace) {
            /* PV_NPlus1(:, xNoPressure, :, 1) =
                    PV_NPlus1(:, xNoPressure, :, 1) .* isNoPressureCell; */
            PV_NPlus1(span, xNoPressure, span, 0) =
                PV_NPlus1(span, xNoPressure, span, 0) * isNoPressureCell;
        }

        // STEP3: Calculate Vx_next, Vy_next and Vz_next
        // CxP [del.Px] = [dPx/dx] = Pr_right - Pr_curr
        // CyP [del.Py] = [dPy/dy] = Pr_top   - Pr_curr
        // CzP [del.Pz] = [dPz/dz] = Pr_front - Pr_curr

        CxP = (PV_NPlus1(seq(1, frameY - 2), seq(2, frameX - 1),
                         seq(1, frameZ - 2), 0) -
               PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 0)) /
              dx;

        CyP = (PV_NPlus1(seq(0, frameY - 3), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 0) -
               PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 0)) /
              dy;

        CzP = (PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(1, frameZ - 2), 0) -
               PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2),
                         seq(2, frameZ - 1), 0)) /
              dz;

        Vx_next(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            (minVxBeta * simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                                      seq(1, frameZ - 2), 1)) -
            (betaVxSqr_dt_invRho * CxP);

        Vy_next(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            (minVyBeta * simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                                      seq(1, frameZ - 2), 2)) -
            (betaVySqr_dt_invRho * CyP);

        Vz_next(seq(0, frameY - 3), seq(0, frameX - 3), seq(0, frameZ - 3)) =
            (minVzBeta * simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                                      seq(1, frameZ - 2), 3)) -
            (betaVzSqr_dt_invRho * CzP);

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  1) = Vx_next(span, span, span);
        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  2) = Vy_next(span, span, span);
        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  3) = Vz_next(span, span, span);

        // STEP4(i) : Inject excitation velocity
        // STEP4(ii): Enforce boundary condition

        // [Note]: Inject excitation velocity as per the source type
        double exeCurrentVal;

        if (sourceModelType == SourceType::VFModel) {
            // TODO
        } else {
            exeCurrentVal = excitationV(T).scalar<float>();
        }

        // Update Vx, Vy and Vz components of the current cell with the
        // excitation velocity
        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  1) +=
            exeCurrentVal * excitation_Vx_weight * maxVxSigmaPrimedt;

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  2) +=
            exeCurrentVal * excitation_Vy_weight * maxVySigmaPrimedt;

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  3) +=
            exeCurrentVal * excitation_Vz_weight * maxVzSigmaPrimedt;

        // Determine velocity near the wall
        vb_alphaX = (vb_alphaX * are_we_not_excitations_Vx.as(f32)) * z_inv;
        vb_alphaY = (vb_alphaY * are_we_not_excitations_Vy.as(f32)) * z_inv;
        vb_alphaZ = (vb_alphaZ * are_we_not_excitations_Vz.as(f32)) * z_inv;

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  1) += maxVxSigmaPrimedt * vb_alphaX;

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  2) += maxVySigmaPrimedt * vb_alphaY;

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  3) += maxVzSigmaPrimedt * vb_alphaZ;

        // Update velocity components
        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  1) /= (minVxBeta + maxVxSigmaPrimedt);

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  2) /= (minVyBeta + maxVySigmaPrimedt);

        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  3) /= (minVzBeta + maxVzSigmaPrimedt);

        // STEP5: Re-store the grid cell type
        PV_NPlus1(seq(1, frameY - 2), seq(1, frameX - 2), seq(1, frameZ - 2),
                  4) = simData.PV_N(seq(1, frameY - 2), seq(1, frameX - 2),
                                    seq(1, frameZ - 2), 4);

        // STEP6: Pass over outer dead cells
        PV_NPlus1(0, span, span, seq(0, 3)) = 0;
        PV_NPlus1(0, span, span, 4) = simData.PV_N(0, span, span, 4);

        PV_NPlus1(frameY - 1, span, span, seq(0, 3)) = 0;
        PV_NPlus1(frameY - 1, span, span, 4) =
            simData.PV_N(frameY - 1, span, span, 4);

        PV_NPlus1(span, 0, span, seq(0, 3)) = 0;
        PV_NPlus1(span, 0, span, 4) = simData.PV_N(span, 0, span, 4);

        PV_NPlus1(span, frameX - 1, span, seq(0, 3)) = 0;
        PV_NPlus1(span, frameX - 1, span, 4) =
            simData.PV_N(span, frameX - 1, span, 4);

        PV_NPlus1(span, span, 0, seq(0, 3)) = 0;
        PV_NPlus1(span, span, 0, 4) = simData.PV_N(span, span, 0, 4);

        PV_NPlus1(span, span, frameZ - 1, seq(0, 3)) = 0;
        PV_NPlus1(span, span, frameZ - 1, 4) =
            simData.PV_N(span, span, frameZ - 1, 4);

        // STEP6: Copy PV_NPlus1to PV_N for the next time frame
        simData.PV_N = PV_NPlus1.copy();

        // Print remaining step numbers
        if (T % 100 == 0) {
            if (plotPressure) {
                // Take slices to visualise the wave propagation

                array slice_xz = simData.PV_N(
                    (uint32_t)std::ceil((frameY - 1) / 2.0), span, span, 0);

                const double minPressure = min<float>(slice_xz);
                const double maxPressure = max<float>(slice_xz);

                slice_xz(simData.PV_N(simData.tubeStart.startY, span, span,
                                      4) == vt::cell_wall) = vis_Boundary;

                const auto slice_xz_image = image::makeWaveImage(
                    moddims(slice_xz, frameX, frameZ), vis_Min, vis_Max);

                image::saveToPng(
                    "images/slice_xz_" + std::to_string(T) + ".png",
                    slice_xz_image);

                fprintf(stderr, "Saved image of pressure at STEP = %d\n", T);
                fprintf(stderr, "  => Min pressure = %f\n", minPressure);
                fprintf(stderr, "  => Max pressure = %f\n", maxPressure);
            }

            fprintf(stderr, "Remaining STEPS = %d\n", STEPS - T);
        }

        // Save audio data as change in pressure
        if (simulationType != SimulationType::OpenSpace) {
            Pr_mouthEnd(T) = PV_NPlus1(simData.listenerInfo.listenerY,
                                       simData.listenerInfo.listenerX,
                                       simData.listenerInfo.listenerZ, 1);

            Pr_glottalEnd(T) = PV_NPlus1(simData.sourceInfo.sourceY,
                                         simData.sourceInfo.sourceX,
                                         simData.sourceInfo.sourceZ, 1);
        }

        // CHECK IF NAN
        bool hasNan = anyTrue<bool>(isNaN(simData.PV_N(span, span, span, 0)));
        if (hasNan) {
            fprintf(stderr, "Solver exploded at step = %d\n", T);
            return;
        }
    }

    const auto listenerAudio =
        audio::generateFromPressure(Pr_mouthEnd, srateBase, srate, srateMul);
    audio::saveToWavFile("listener.wav", listenerAudio, srate);

    const auto sourceAudio =
        audio::generateFromPressure(Pr_glottalEnd, srateBase, srate, srateMul);
    audio::saveToWavFile("source.wav", sourceAudio, srate);
}
#include "audio/generate_audio.hh"
#include "audio/save_audio_file.hh"
#include "constants.hh"
#include "image/make_cell_image.hh"
#include "image/save_to_png.hh"
#include "routines/circular_tube_generation.hh"
#include "routines/elliptical_tube_generation.hh"
#include "routines/impulse_signal.hh"

using namespace vt;

void talkingTube() {
    const int              srateBase = 44100;
    const uint32_t         srateMul = 10;  // Sample rate multiplier
    const double           dur = 50e-3;    // simulation duration (50 ms)
    const SimulationType   simulationType = SimulationType::VowelSound;
    const Vowel            vowel = Vowel::i;
    const CrossSectionType crossSectionType = CrossSectionType::Elliptical;
    const JunctionType     junctionType = JunctionType::Centric;
    const SourceType       sourceModelType = SourceType::Impulse;
    const bool             pmlSwitch = false;
    const bool             baffleSwitchFlag = false;
    const MouthTermination rad = MouthTermination::DirichletBoundary;
    const bool             boundaryInterpolation = false;

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
        curTime = timeIndex * dt;
        timeIndex++;
    }
    const uint32_t STEPS = times.size();
    auto           t = zeros<ArrayXd>(STEPS);
    std::copy(times.begin(), times.end(), std::next(t.begin()));

    // DEFINE BETA AND SIGMAPRIME PARAMETERS

    const uint32_t vis_Boundary = 2000;
    auto           sigmadt = zeros<ArrayXd>(pmlLayer);

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
    auto typeValues = zeros<Array2Xd>(2, vt::cell_numTypes);
    typeValues(seq(1, 2), vt::cell_wall + 1) << 0, dt;        // VT walls
    typeValues(seq(1, 2), vt::cell_air + 1) << 1, 0;          // air
    typeValues(seq(1, 2), vt::cell_noPressure + 1) << 1, 0;   // air
    typeValues(seq(1, 2), vt::cell_excitation + 1) << 0, dt;  // excitation
    typeValues(seq(1, 2), vt::cell_head + 1) << 0, dt;        // head cells

    // Define beta and sigmaPrimedt for PML layers
    // For PML layers beta = 1,
    // sigmaPrimedt = [1-beta+sigma]*dt = [1-1+sigma]*dt = sigmadt
    for (uint32_t pmlCounter = 0; pmlCounter <= pmlLayer - 1; ++pmlCounter) {
        sigmadt(pmlCounter + 1) = (pmlCounter / (pmlLayer - 1)) * maxSigmadt;
        typeValues(seq(1, 2), vt::cell_pml0 + 1 + pmlCounter) << 1,
            sigmadt(pmlCounter + 1);
    }
    typeValues(seq(1, 2), vt::cell_dead + 1) << 0, 1000000;

    // DEFINE SIMULATION TYPE

    // [Note]: I have created two separate functions for the circular and
    // elliptical tube generation. However, that's not required since circle is
    // a special case of ellipse where the ratio between the semiMajorAxis and
    // semiMinorAxis is equal. [semiMajorAxis:semiMinorAxis = 1:1]

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

    auto Pr_mouthEnd = zeros<ArrayXd>(STEPS);
    auto Pr_glottalEnd = zeros<ArrayXd>(STEPS);

    // DETERMINE WALL IMPEDANCE

    const double mu3D = 0.005;  // Boundary Admittance Coefficient
    const double alpha3D =
        1 / (0.5 + 0.25 * (mu3D + (1 / mu3D)));  // Sound Absorption Coefficient
    const double z_inv =
        1 / (rhoC * ((1 + sqrt(1 - alpha3D)) / (1 - sqrt(1 - alpha3D))));

    // SOURCE MODEL

    auto excitationV = zeros<ArrayXd>(STEPS);

    switch (sourceModelType) {
        case SourceType::Sine: {
            // Sine wave source model
            const double excitationF = 440;  // 440 Hz
            const double srcAmplitude = 25;
            for (uint32_t i = 1; i <= STEPS; ++i) {
                excitationV(i) =
                    srcAmplitude * sin(2 * pi * excitationF * dt * (i - 1));
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
    auto srcDirection = zeros<ArrayXd>(6);
    // For forward direction
    srcDirection(4) = 1;

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
    for (uint32_t x = 1; x <= frameX; ++x) {
        for (uint32_t z = 1; z <= frameZ; ++z) {
            simData.PV_N(1, x, z, 5) = vt::cell_dead;
            simData.PV_N(frameY, x, z, 5) = vt::cell_dead;
        }
    }

    // ---Adding dead cells to the back and front surfaces---
    for (uint32_t y = 1; y <= frameY; ++y) {
        for (uint32_t x = 1; x <= frameX; ++x) {
            simData.PV_N(y, x, 1, 5) = vt::cell_dead;
            simData.PV_N(y, x, frameZ, 5) = vt::cell_dead;
        }
    }

    // ---Adding dead cells to the left and right surfaces---
    for (uint32_t y = 1; y <= frameY; ++y) {
        for (uint32_t z = 1; z <= frameZ; ++z) {
            simData.PV_N(y, 1, z, 5) = vt::cell_dead;
            simData.PV_N(y, frameX, z, 5) = vt::cell_dead;
        }
    }

    if (pmlSwitch) {
        // Define PML layers starting from the outer-most layer
        uint32_t pmlType = (uint32_t)vt::cell_pml5;

        const uint32_t xShift = 1;
        uint32_t       xStart = 2;
        uint32_t       xEnd = frameX - 1;

        const uint32_t yShift = 1;
        uint32_t       yStart = 2;
        uint32_t       yEnd = frameY - 1;

        const uint32_t zShift = 1;
        uint32_t       zStart = 2;
        uint32_t       zEnd = frameZ - 1;

        for (uint32_t pmlCount = 1; pmlCount <= pmlLayer; ++pmlCount) {
            // ---Adding pml layers to the top and bottom surfaces---
            for (uint32_t x = xStart; x <= xEnd; ++x) {
                for (uint32_t z = zStart; z <= zEnd; ++z) {
                    simData.PV_N(yShift + pmlCount, x, z, 5) = pmlType;
                    simData.PV_N(frameY - pmlCount, x, z, 5) = pmlType;
                }
            }

            // ---Adding pml layers to the back and front surfaces---
            for (uint32_t y = yStart; y <= yEnd; ++y) {
                for (uint32_t x = xStart; x <= xEnd; ++x) {
                    simData.PV_N(y, x, zShift + pmlCount, 5) = pmlType;
                    simData.PV_N(y, x, frameZ - pmlCount, 5) = pmlType;
                }
            }

            // ---Adding pml layers to the left and right surfaces---
            for (uint32_t y = yStart; y <= yEnd; ++y) {
                for (uint32_t z = zStart; z <= zEnd; ++z) {
                    simData.PV_N(y, xShift + pmlCount, z, 5) = pmlType;
                    simData.PV_N(y, frameX - pmlCount, z, 5) = pmlType;
                }
            }

            // Update the pmlLayer type.
            // Update the start and end position
            pmlType--;

            xStart++;
            xEnd--;

            yStart++;
            yEnd--;

            zStart++;
            zEnd--;
        }
    }

    // MODEL VISUALISATION

    // Visualize the PML layers and excitation cell by taking a slice of the
    // domain
    Tensor3<GridCellType> buildFrame =
        simData.PV_N.chip(5, 3).cast<GridCellType>();
    // Uncomment to visualise source position for the open space simulation
    // buildFrame(midY, midX, midZ) = -1;
    buildFrame(simData.listenerInfo.listenerY, simData.listenerInfo.listenerX,
               simData.listenerInfo.listenerZ) = GridCellType::cell_head;

    // Slice along XZ-plane
    const auto frame_xz = buildFrame.chip((uint32_t)std::ceil(frameY / 2.0), 0);

    // Slice along YZ-plane
    const auto frame_yz = buildFrame.chip((uint32_t)std::ceil(frameX / 2.0), 1);

    // Slice along XY-plane
    const auto frame_xy = buildFrame.chip((uint32_t)std::ceil(frameZ / 2.0), 2);

    // Save images
    const auto frame_xz_image = image::makeCellImage(frame_xz);
    const auto frame_yz_image = image::makeCellImage(frame_yz);
    const auto frame_xy_image = image::makeCellImage(frame_xy);
    image::saveToPng("frame_xz.png", frame_xz_image);
    image::saveToPng("frame_yz.png", frame_yz_image);
    image::saveToPng("frame_xy.png", frame_xy_image);

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
    auto cellTypes = zeros<Array4<uint32_t>>(4);
    auto typeIndex = zeros<Array4<uint32_t>>(4);
    auto beta = zeros<Array4<double>>(4);
    auto sigmaPrime_dt = zeros<Array4<double>>(4);

    // To store is_excitation and excitation_weight
    const uint32_t numGridCellsForComputation = height * width * depth;
    auto is_excitation = zeros<Array4X<bool>>(4, numGridCellsForComputation);
    auto excitation_weight =
        zeros<Array3X<double>>(3, numGridCellsForComputation);
    auto N_out = zeros<ArrayX<double>>(numGridCellsForComputation);
    auto N_in = zeros<ArrayX<double>>(numGridCellsForComputation);

    // Create arrays to store minBeta and maxSigmaPrime_dt
    // Note => We consider minBeta, because beta parameter is defined at the
    // centre of the cell. And the velocity components are defined either on
    // edges [2D grid cell]/ side surfaces [3D grid cell]. Therefore, we
    // tie-break between these grid cells considering minBeta value.

    auto minVxBeta = zeros<Tensor3<double>>(height, width, depth);
    auto minVyBeta = zeros<Tensor3<double>>(height, width, depth);
    auto minVzBeta = zeros<Tensor3<double>>(height, width, depth);

    auto maxVxSigmaPrimedt = zeros<Tensor3<double>>(height, width, depth);
    auto maxVySigmaPrimedt = zeros<Tensor3<double>>(height, width, depth);
    auto maxVzSigmaPrimedt = zeros<Tensor3<double>>(height, width, depth);

    // Store the excitation_weight and are_we_not_excitations in matrix format
    // to compute Vx, Vy, Vz
    auto excitation_Vx_weight = zeros<Tensor3<double>>(height, width, depth);
    auto excitation_Vy_weight = zeros<Tensor3<double>>(height, width, depth);
    auto excitation_Vz_weight = zeros<Tensor3<double>>(height, width, depth);

    auto are_we_not_excitations_Vx = zeros<Tensor3<bool>>(height, width, depth);
    auto are_we_not_excitations_Vy = zeros<Tensor3<bool>>(height, width, depth);
    auto are_we_not_excitations_Vz = zeros<Tensor3<bool>>(height, width, depth);

    // Store the xor_val in matrix format
    auto xor_val1 = zeros<Tensor3<bool>>(height, width, depth);
    auto xor_val2 = zeros<Tensor3<bool>>(height, width, depth);
    auto xor_val3 = zeros<Tensor3<bool>>(height, width, depth);
    auto xor_val4 = zeros<Tensor3<bool>>(height, width, depth);
    auto xor_val5 = zeros<Tensor3<bool>>(height, width, depth);
    auto xor_val6 = zeros<Tensor3<bool>>(height, width, depth);

    // Store the N_out and N_in in matrix format
    auto N_out_mat = zeros<Tensor3<double>>(height, width, depth);
    auto N_in_mat = zeros<Tensor3<double>>(height, width, depth);

    // Note => We are reshaping sigmaPrimedt in matrix format as we need this
    // while calculating pressure [Check the denominator of the discretized
    // pressure equation]
    auto sigmaPrimedt = zeros<Tensor3<double>>(height, width, depth);

    // Note => Here for each grid cell, we store cell type, beta and
    // siggmaPrine_dt values in arrays for its neighbor cells and the current
    // cell. In case of For example, in case of -
    // 2D/2.5D => [current, right_current, top_current]
    // 3D => [current, right_current, up_current, front_current]

    // Set a counter to save data for is_excitation, excitation_weight, xor,
    // N_out, N_in
    uint32_t counter = 1;

    for (uint32_t height_idx = 2; height_idx <= frameY - 1; ++height_idx) {
        for (uint32_t width_idx = 2; width_idx <= frameX - 1; ++width_idx) {
            for (uint32_t depth_idx = 2; depth_idx <= frameZ - 1; ++depth_idx) {
                // Find the cellTypes = [current, right_current, up_current,
                // front_current]
                auto cellTypes = zeros<Array4<uint32_t>>(4);
                cellTypes(1) =
                    simData.PV_N(height_idx, width_idx, depth_idx, 5);
                cellTypes(2) =
                    simData.PV_N(height_idx, width_idx + 1, depth_idx, 5);
                cellTypes(3) =
                    simData.PV_N(height_idx - 1, width_idx, depth_idx, 5);
                cellTypes(4) =
                    simData.PV_N(height_idx, width_idx, depth_idx + 1, 5);

                // For typeIndex add 1 to cellTypes
                typeIndex = cellTypes + 1;

                // Store beta values in beta array
                for (uint32_t i = 1; i <= 4; ++i) {
                    beta(i) = typeValues(1, typeIndex(i));
                }

                // Calculate minBeta
                const double min_beta_Vx = std::min(
                    beta(1), beta(2));  // minBeta(current, right_current)
                const double min_beta_Vy = std::min(
                    beta(1), beta(3));  // minBeta(current, top_current)
                const double min_beta_Vz = std::min(
                    beta(1), beta(4));  // minBeta(current, front_current)

                minVxBeta(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    min_beta_Vx;
                minVyBeta(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    min_beta_Vy;
                minVzBeta(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    min_beta_Vz;

                // Store sigmaPrime*dt values in sigmaPrime_dt
                for (uint32_t i = 1; i <= 4; ++i) {
                    sigmaPrime_dt(i) = typeValues(2, typeIndex(i));
                }

                // Store sigmaPrime_dt of the current cell only
                sigmaPrimedt(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    sigmaPrime_dt(1);

                // Calculate maxSigmaPrime_dt
                const double max_sigmaPrimedt_Vx =
                    std::max(sigmaPrime_dt(1), sigmaPrime_dt(2));
                const double max_sigmaPrimedt_Vy =
                    std::max(sigmaPrime_dt(1), sigmaPrime_dt(3));
                const double max_sigmaPrimedt_Vz =
                    std::max(sigmaPrime_dt(1), sigmaPrime_dt(4));

                maxVxSigmaPrimedt(height_idx - 1, width_idx - 1,
                                  depth_idx - 1) = max_sigmaPrimedt_Vx;
                maxVySigmaPrimedt(height_idx - 1, width_idx - 1,
                                  depth_idx - 1) = max_sigmaPrimedt_Vy;
                maxVzSigmaPrimedt(height_idx - 1, width_idx - 1,
                                  depth_idx - 1) = max_sigmaPrimedt_Vz;

                // Check whether the current cell is an excitation cell or not
                for (uint32_t i = 1; i <= 4; ++i) {
                    is_excitation(i, counter) =
                        (cellTypes(i) == (uint32_t)vt::cell_excitation);
                }

                // Verify if both the current cell and neighbouring cells are
                // not excitation cells
                are_we_not_excitations_Vx(height_idx - 1, width_idx - 1,
                                          depth_idx - 1) =
                    (!is_excitation(1, counter) && !is_excitation(2, counter));

                are_we_not_excitations_Vy(height_idx - 1, width_idx - 1,
                                          depth_idx - 1) =
                    (!is_excitation(1, counter) && !is_excitation(3, counter));

                are_we_not_excitations_Vz(height_idx - 1, width_idx - 1,
                                          depth_idx - 1) =
                    (!is_excitation(1, counter) && !is_excitation(4, counter));

                // Determine excitation weight for the current cell
                // Find excitation velocity going out of the cell and coming
                // back into the cells depending upon the source direction.
                // Add them all to find the net velocity [velocity components]
                // for the current cell.

                const auto excitationWeightForward =
                    is_excitation(1, counter) * (srcDirection(seq(4, 6)));

                const auto excitationWeightBackward =
                    is_excitation(seq(2, 4), counter).cast<double>() *
                    (srcDirection(seq(1, 3)));

                excitation_weight(seq(1, 3), counter) =
                    excitationWeightForward + excitationWeightBackward;

                // Store the excitation_weight in matrix format
                excitation_Vx_weight(height_idx - 1, width_idx - 1,
                                     depth_idx - 1) =
                    excitation_weight(1, counter);

                excitation_Vy_weight(height_idx - 1, width_idx - 1,
                                     depth_idx - 1) =
                    excitation_weight(2, counter);

                excitation_Vz_weight(height_idx - 1, width_idx - 1,
                                     depth_idx - 1) =
                    excitation_weight(3, counter);

                // Check the adjacent cells to determine the velocity direction
                xor_val1(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(2) * (1 - beta(1));
                xor_val2(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(1) * (1 - beta(2));

                xor_val3(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(3) * (1 - beta(1));
                xor_val4(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(1) * (1 - beta(3));

                xor_val5(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(4) * (1 - beta(1));
                xor_val6(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    beta(1) * (1 - beta(4));

                // Find the boundary condition from the lookup table and
                // select the corresponding unit vector
                uint32_t lookup_idx = (uint32_t)(-1);

                for (uint32_t i = 0; i < lookUp_table_count; ++i) {
                    if (lookUp_table[i][0] == beta(2) &&
                        lookUp_table[i][1] == beta(3) &&
                        lookUp_table[i][2] == beta(4)) {
                        lookup_idx = i;
                        break;
                    }
                }

                if (lookup_idx == (uint32_t)(-1)) {
                    throw std::logic_error("invalid state");
                }

                N_out(counter) = air_normalV_component[lookup_idx] * beta(1);
                N_in(counter) =
                    wall_normalV_component[lookup_idx] * (1 - beta(1));

                N_out_mat(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    N_out(counter);
                N_in_mat(height_idx - 1, width_idx - 1, depth_idx - 1) =
                    N_in(counter);

                // Increment the counter
                counter++;
            }
        }
    }

    const auto betaVxSqr = minVxBeta * minVxBeta;
    const auto betaVySqr = minVyBeta * minVyBeta;
    const auto betaVzSqr = minVzBeta * minVzBeta;

    const Tensor3<double> betaVxSqr_dt_invRho = (betaVxSqr * dt) / rho;
    const Tensor3<double> betaVySqr_dt_invRho = (betaVySqr * dt) / rho;
    const Tensor3<double> betaVzSqr_dt_invRho = (betaVzSqr * dt) / rho;

    const auto rho_sqrC_dt_invds = (kappa * dt) / ds;
    const auto rho_sqrC_dt = kappa * dt;

    // INITIALISE ACOUSTIC PARAMETERS AND COEFFICIENTS
    auto PV_NPlus1 = zeros<Tensor4<double>>(frameY, frameX, frameZ, 5);

    auto Ug_array = zeros<ArrayXd>(STEPS);

    auto CxVx = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto CyVy = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto CzVz = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);

    auto CxP = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto CyP = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto CzP = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);

    auto Pr_next = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto Vx_next = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto Vy_next = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto Vz_next = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);

    auto vb_alphaX = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto vb_alphaY = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);
    auto vb_alphaZ = zeros<Tensor3<double>>(frameY - 2, frameX - 2, frameZ - 2);

    // Introduce Dirichlet boundary condition
    uint32_t      xNoPressure;
    Tensor2<bool> isNoPressureCell;

    if (rad == MouthTermination::DirichletBoundary &&
        simulationType != SimulationType::OpenSpace) {
        xNoPressure = simData.tubeStart.startX + simData.totalTubeLengthInCells;
        // gridNoPressurePlane = PV_N(:, xNoPressure, :, 5)
        const auto gridNoPressurePlane =
            simData.PV_N.chip(5, 3).chip(xNoPressure, 1);
        isNoPressureCell =
            (gridNoPressurePlane == GridCellType::cell_noPressure);
    }

    for (uint32_t T = 1; T <= STEPS; ++T) {
        // STEP1: Calculate [del.V] = [dVx/dx + dVy/dy + dVz/dz]
        // CxVx = dVx, where Vx = velocity along x direction = Vx_curr - Vx_left
        // CyVy = dVy, where Vy = velocity along y direction = Vy_curr - Vy_down
        // CzVz = dVz, where Vz = velocity along z direction = Vz_curr - Vz_back

        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    CxVx(y, x, z) = simData.PV_N(y + 1, x + 1, z + 1, 2) -
                                    simData.PV_N(y + 1, x, z + 1, 2);

                    CyVy(y, x, z) = simData.PV_N(y + 1, x + 1, z + 1, 3) -
                                    simData.PV_N(y + 2, x + 1, z + 1, 3);

                    CzVz(y, x, z) = simData.PV_N(y + 1, x + 1, z + 1, 4) -
                                    simData.PV_N(y + 1, x + 1, z, 4);
                }
            }
        }

        // STEP2: Calculate Pr_next

        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    Pr_next(y, x, z) =
                        (simData.PV_N(y + 1, x + 1, z + 1, 1) -
                         (rho_sqrC_dt_invds *
                          (CxVx(y, x, z) + CyVy(y, x, z) + CzVz(y, x, z)))) /
                        (1 + sigmaPrimedt(y, x, z));

                    PV_NPlus1(y + 1, x + 1, z + 1, 1) = Pr_next(y, x, z);
                }
            }
        }

        // Make pressure values of no_pressure cells as zeros
        if (rad == MouthTermination::DirichletBoundary &&
            simulationType != SimulationType::OpenSpace) {
            /* PV_NPlus1(:, xNoPressure, :, 1) =
                    PV_NPlus1(:, xNoPressure, :, 1) .* isNoPressureCell; */
            PV_NPlus1.chip(1, 3).chip(xNoPressure, 1) *=
                isNoPressureCell.cast<double>();
        }

        // STEP3: Calculate Vx_next, Vy_next and Vz_next
        // CxP [del.Px] = [dPx/dx] = Pr_right - Pr_curr
        // CyP [del.Py] = [dPy/dy] = Pr_top   - Pr_curr
        // CzP [del.Pz] = [dPz/dz] = Pr_front - Pr_curr

        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    CxP(y, x, z) = (PV_NPlus1(y + 1, x + 2, z + 2, 1) -
                                    PV_NPlus1(y + 1, x + 1, z + 1, 1)) /
                                   dx;

                    CyP(y, x, z) = (PV_NPlus1(y, x + 1, z + 1, 1) -
                                    PV_NPlus1(y + 1, x + 1, z + 1, 1)) /
                                   dy;

                    CzP(y, x, z) = (PV_NPlus1(y + 1, x + 1, z + 2, 1) -
                                    PV_NPlus1(y + 1, x + 1, z + 1, 1)) /
                                   dz;

                    Vx_next(y, x, z) =
                        (minVxBeta(y, x, z) *
                         simData.PV_N(y + 1, x + 1, z + 1, 2)) -
                        (betaVxSqr_dt_invRho(y, x, z) * CxP(y, x, z));

                    Vy_next(y, x, z) =
                        (minVyBeta(y, x, z) *
                         simData.PV_N(y + 1, x + 1, z + 1, 3)) -
                        (betaVySqr_dt_invRho(y, x, z) * CyP(y, x, z));

                    Vz_next(y, x, z) =
                        (minVzBeta(y, x, z) *
                         simData.PV_N(y + 1, x + 1, z + 1, 4)) -
                        (betaVzSqr_dt_invRho(y, x, z) * CzP(y, x, z));

                    PV_NPlus1(y + 1, x + 1, z + 1, 2) = Vx_next(y, x, z);
                    PV_NPlus1(y + 1, x + 1, z + 1, 3) = Vy_next(y, x, z);
                    PV_NPlus1(y + 1, x + 1, z + 1, 4) = Vz_next(y, x, z);
                }
            }
        }

        // STEP4(i) : Inject excitation velocity
        // STEP4(ii): Enforce boundary condition

        // [Note]: Inject excitation velocity as per the source type
        double exeCurrentVal;

        if (sourceModelType == SourceType::VFModel) {
            // TODO
        } else {
            exeCurrentVal = excitationV(T);
        }

        // Update Vx, Vy and Vz components of the current cell with the
        // excitation velocity
        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    PV_NPlus1(y + 1, x + 1, z + 1, 2) +=
                        exeCurrentVal * excitation_Vx_weight(y, x, z) *
                        maxVxSigmaPrimedt(y, x, z);

                    PV_NPlus1(y + 1, x + 1, z + 1, 3) +=
                        exeCurrentVal * excitation_Vy_weight(y, x, z) *
                        maxVySigmaPrimedt(y, x, z);

                    PV_NPlus1(y + 1, x + 1, z + 1, 4) +=
                        exeCurrentVal * excitation_Vz_weight(y, x, z) *
                        maxVzSigmaPrimedt(y, x, z);
                }
            }
        }

        // Determine velocity near the wall
        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    vb_alphaX(y, x, z) =
                        xor_val2(y, x, z) * PV_NPlus1(y + 1, x + 1, z + 1, 1) *
                            N_out_mat(y, x, z) -
                        xor_val1(y, x, z) * PV_NPlus1(y + 1, x + 2, z + 1, 1) *
                            N_in_mat(y, x, z);

                    vb_alphaY(y, x, z) =
                        xor_val4(y, x, z) * PV_NPlus1(y + 1, x + 1, z + 1, 1) *
                            N_out_mat(y, x, z) -
                        xor_val3(y, x, z) * PV_NPlus1(y, x + 1, z + 1, 1) *
                            N_in_mat(y, x, z);

                    vb_alphaZ(y, x, z) =
                        xor_val6(y, x, z) * PV_NPlus1(y + 1, x + 1, z + 1, 1) *
                            N_out_mat(y, x, z) -
                        xor_val5(y, x, z) * PV_NPlus1(y + 1, x + 1, z + 2, 1) *
                            N_in_mat(y, x, z);
                }
            }
        }

        // Apply tube boundary condition/wall losses
        vb_alphaX *= are_we_not_excitations_Vx.cast<double>() * z_inv;
        vb_alphaY *= are_we_not_excitations_Vy.cast<double>() * z_inv;
        vb_alphaZ *= are_we_not_excitations_Vz.cast<double>() * z_inv;

        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    PV_NPlus1(y + 1, x + 1, z + 1, 2) +=
                        maxVxSigmaPrimedt(y, x, z) * vb_alphaX(y, x, z);

                    PV_NPlus1(y + 1, x + 1, z + 1, 3) +=
                        maxVySigmaPrimedt(y, x, z) * vb_alphaY(y, x, z);

                    PV_NPlus1(y + 1, x + 1, z + 1, 4) +=
                        maxVzSigmaPrimedt(y, x, z) * vb_alphaZ(y, x, z);
                }
            }
        }

        // Update velocity components
        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    PV_NPlus1(y + 1, x + 1, z + 1, 2) /=
                        (minVxBeta(y, x, z) + maxVxSigmaPrimedt(y, x, z));

                    PV_NPlus1(y + 1, x + 1, z + 1, 3) /=
                        (minVyBeta(y, x, z) + maxVySigmaPrimedt(y, x, z));

                    PV_NPlus1(y + 1, x + 1, z + 1, 4) /=
                        (minVzBeta(y, x, z) + maxVzSigmaPrimedt(y, x, z));
                }
            }
        }

        // STEP5: Re-store the grid cell type
        for (uint32_t y = 1; y <= frameY - 2; ++y) {
            for (uint32_t x = 1; x <= frameX - 2; ++x) {
                for (uint32_t z = 1; z <= frameZ - 2; ++z) {
                    PV_NPlus1(y + 1, x + 1, z + 1, 5) =
                        simData.PV_N(y + 1, x + 1, z + 1, 5);
                }
            }
        }

        // STEP6: Pass over outer dead cells
        for (uint32_t x = 1; x <= frameX; ++x) {
            for (uint32_t z = 1; z <= frameZ; ++z) {
                for (uint32_t k = 1; k <= 4; ++k) {
                    PV_NPlus1(1, x, z, k) = 0;
                    PV_NPlus1(frameY, x, z, k) = 0;
                }
                PV_NPlus1(1, x, z, 5) = simData.PV_N(1, x, z, 5);
                PV_NPlus1(frameY, x, z, 5) = simData.PV_N(frameY, x, z, 5);
            }
        }

        for (uint32_t y = 1; y <= frameY; ++y) {
            for (uint32_t z = 1; z <= frameZ; ++z) {
                for (uint32_t k = 1; k <= 4; ++k) {
                    PV_NPlus1(y, 1, z, k) = 0;
                    PV_NPlus1(y, frameX, z, k) = 0;
                }
                PV_NPlus1(y, 1, z, 5) = simData.PV_N(y, 1, z, 5);
                PV_NPlus1(y, frameX, z, 5) = simData.PV_N(y, frameX, z, 5);
            }
        }

        for (uint32_t y = 1; y <= frameY; ++y) {
            for (uint32_t x = 1; x <= frameX; ++x) {
                for (uint32_t k = 1; k <= 4; ++k) {
                    PV_NPlus1(y, x, 1, k) = 0;
                    PV_NPlus1(y, x, frameZ, k) = 0;
                }
                PV_NPlus1(y, x, 1, 5) = simData.PV_N(y, x, 1, 5);
                PV_NPlus1(y, x, frameZ, 5) = simData.PV_N(y, x, frameZ, 5);
            }
        }

        // STEP6: Copy PV_NPlus1to PV_N for the next time frame
        simData.PV_N = PV_NPlus1;

        // Print remaining step numbers
        if (T % 1000 == 0) {
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
        bool hasNan = false;

        for (uint32_t z = 1; z <= frameZ; ++z) {
            for (uint32_t x = 1; x <= frameX; ++x) {
                const double test_xz =
                    simData.PV_N((uint32_t)std::ceil(frameY / 2.0), x, z, 1);
                if (std::isnan(test_xz)) {
                    hasNan = true;
                    break;
                }
            }
            if (hasNan) break;

            for (uint32_t y = 1; y <= frameY; ++y) {
                const double test_yz =
                    simData.PV_N(y, (uint32_t)std::ceil(frameX / 2.0), z, 1);
                if (std::isnan(test_yz)) {
                    hasNan = true;
                    break;
                }
            }
            if (hasNan) break;
        }

        if (hasNan) {
            fprintf(stderr, "Solver exploded at step = %d\n", T);
            return;
        }
    }

    const auto listenerAudio =
        audio::generateFromPressure(Pr_mouthEnd, srateBase, srate, srateMul);
    audio::saveToWavFile("listener.wav", listenerAudio, srateBase);

    const auto sourceAudio =
        audio::generateFromPressure(Pr_mouthEnd, srateBase, srate, srateMul);
    audio::saveToWavFile("source.wav", sourceAudio, srateBase);
}
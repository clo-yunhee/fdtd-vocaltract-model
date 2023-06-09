#include "area_function.hh"

static vt::TubeProperties s_regularTube(
    {
        0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
        0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
        0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
        0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97,
    },      // in cm^2
    0.1927  // in m
);

static vt::TubeProperties s_vowel_a(
    {
        0.56, 0.62, 0.66, 0.78, 0.97, 1.16, 1.12, 0.82, 0.55, 0.45, 0.37,
        0.29, 0.21, 0.15, 0.16, 0.25, 0.34, 0.43, 0.54, 0.61, 0.67, 0.98,
        1.76, 2.75, 3.52, 4.08, 4.74, 5.61, 6.60, 7.61, 8.48, 9.06, 9.29,
        9.26, 9.06, 8.64, 7.91, 6.98, 6.02, 5.13, 4.55, 4.52, 4.71, 4.72,
    },      // in cm^2
    0.1709  // in m
);

static vt::TubeProperties s_vowel_i(
    {
        0.51, 0.59, 0.62, 0.72, 1.24, 2.30, 3.30, 3.59, 3.22, 2.86, 3.00,
        3.61, 4.39, 4.95, 5.17, 5.16, 5.18, 5.26, 5.20, 5.02, 4.71, 4.13,
        3.43, 2.83, 2.32, 1.83, 1.46, 1.23, 1.08, 0.94, 0.80, 0.67, 0.55,
        0.46, 0.40, 0.36, 0.35, 0.35, 0.38, 0.51, 0.74, 0.92, 0.96, 0.91,
    },      // in cm^2
    0.1690  // in m
);

static vt::TubeProperties s_vowel_u(
    {
        0.54, 0.61, 0.66, 0.75, 1.13, 1.99, 2.83, 2.90, 2.52, 2.40, 2.83,
        3.56, 3.99, 3.89, 3.50, 3.04, 2.64, 2.44, 2.31, 2.07, 1.80, 1.52,
        1.14, 0.74, 0.42, 0.22, 0.14, 0.20, 0.47, 0.89, 1.15, 1.42, 2.17,
        3.04, 3.69, 4.70, 5.74, 5.41, 3.82, 2.34, 1.35, 0.65, 0.29, 0.16,
    },      // in cm^2
    0.1959  // in m
);

static vt::TubeProperties s_vowel_e(
    {
        0.29, 0.26, 0.30, 0.40, 0.55, 0.74, 0.99, 1.09, 0.90, 0.69, 0.77,
        1.31, 2.13, 2.74, 3.03, 3.23, 3.33, 3.27, 3.09, 2.84, 2.66, 2.46,
        2.14, 1.79, 1.44, 1.17, 1.00, 0.88, 0.80, 0.81, 0.85, 0.84, 0.86,
        0.96, 1.18, 1.35, 1.48, 1.62, 1.49, 1.29, 1.24, 1.17, 1.04, 0.95,
    },      // in cm^2
    0.1698  // in m
);

static vt::TubeProperties s_vowel_o(
    {
        0.38, 0.45, 0.57, 0.77, 1.31, 1.92, 1.74, 1.11, 0.75, 0.59, 0.57,
        0.68, 0.73, 0.67, 0.58, 0.49, 0.44, 0.42, 0.49, 0.53, 0.38, 0.30,
        0.45, 0.61, 0.71, 0.79, 0.86, 1.01, 1.41, 2.09, 3.00, 4.10, 5.16,
        6.22, 7.34, 8.15, 8.61, 8.37, 6.76, 4.37, 2.30, 1.06, 0.58, 0.47,
    },      // in cm^2
    0.1833  // in m
);

static vt::TubeProperties s_vowel_I(
    {
        0.28, 0.21, 0.21, 0.30, 0.47, 0.71, 1.12, 1.48, 1.35, 1.05, 0.92,
        0.92, 1.19, 1.94, 2.83, 3.31, 3.48, 3.60, 3.64, 3.49, 3.20, 2.90,
        2.59, 2.21, 1.87, 1.54, 1.20, 0.92, 0.74, 0.59, 0.52, 0.54, 0.59,
        0.65, 0.71, 0.67, 0.61, 0.57, 0.50, 0.48, 0.54, 0.73, 0.93, 0.82,
    },      // in cm^2
    0.1655  // in m
);

vt::TubeProperties vt::areaFunction(SimulationType simulationType,
                                    JunctionType junctionType, Vowel vowel) {
    vt::TubeProperties tube{s_regularTube};

    if (simulationType != SimulationType::RegularTube) {
        switch (vowel) {
            case Vowel::a:
                tube = s_vowel_a;
                break;
            case Vowel::i:
                tube = s_vowel_i;
                break;
            case Vowel::u:
                tube = s_vowel_u;
                break;
            case Vowel::e:
                tube = s_vowel_e;
                break;
            case Vowel::o:
                tube = s_vowel_o;
                break;
            case Vowel::I:
                tube = s_vowel_I;
                break;
            default:
                throw std::invalid_argument(
                    "Vowel can't be None if simulationType isn't RegularType");
        }
    }

    if (junctionType != JunctionType::Centric) {
        switch (vowel) {
            case Vowel::a:
                tube.segmentLength *= 0.932;
                break;
            case Vowel::i:
                tube.segmentLength *= 0.932;
                break;
            case Vowel::u:
                tube.segmentLength *= 0.89;
                break;
            case Vowel::e:
                tube.segmentLength *= 0.952;
                break;
            case Vowel::o:
                tube.segmentLength *= 0.88;
                break;
            default:
                break;
        }
    }

    return tube;
}
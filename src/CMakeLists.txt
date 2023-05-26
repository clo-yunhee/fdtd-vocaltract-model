set(_target VocalTractModel)

add_executable(${_target}
    audio/butterworth.cc
    audio/butterworth.hh
    audio/generate_audio.cc
    audio/generate_audio.hh
    audio/save_audio_file.cc
    audio/save_audio_file.hh
    audio/sos_filter.cc
    audio/sos_filter.hh
    image/image.hh
    image/make_cell_image.cc
    image/make_cell_image.hh
    image/make_wave_image.cc
    image/make_wave_image.hh
    image/save_to_png.cc
    image/save_to_png.hh
    model/LF_precomputed_Rd_double.inc.h
    model/LF.cpp
    model/LF.h
    routines/area_function.cc
    routines/area_function.hh
    routines/boundary_interpolation.cc
    routines/boundary_interpolation.hh
    routines/centric_circular_contour_generation.cc
    routines/centric_circular_contour_generation.hh
    routines/centric_elliptical_contour_generation.cc
    routines/centric_elliptical_contour_generation.hh
    routines/circular_tube_generation.cc
    routines/circular_tube_generation.hh
    routines/eccentric_circular_contour_generation.cc
    routines/eccentric_circular_contour_generation.hh
    routines/eccentric_elliptical_contour_generation.cc
    routines/eccentric_elliptical_contour_generation.hh
    routines/elliptical_tube_generation.cc
    routines/elliptical_tube_generation.hh
    routines/find_cell_types.cc
    routines/find_cell_types.hh
    routines/impulse_signal.cc
    routines/impulse_signal.hh
    routines/tube_segment_connector.cc
    routines/tube_segment_connector.hh
    main.cc
    talking_tube.cc
    types.hh
)

target_compile_features(${_target} PUBLIC cxx_std_20)
set_target_properties(${_target} PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(${_target} PRIVATE
    ArrayFire::af
    SndFile::sndfile
    PNG::PNG
    Boost::headers
    colormap::colormap
)

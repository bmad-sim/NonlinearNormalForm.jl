cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set(INC_DIRS 
)

set (SRC_DIRS
 .
)

set (EXE_SPECS
  cmake_files/cmake.z_deriv
  cmake_files/cmake.z_sub_i
  cmake_files/cmake.z_cut
  cmake_files/cmake.z_sub_j
  cmake_files/cmake.z_var
  cmake_files/cmake.z_mono
  cmake_files/cmake.z_par
  cmake_files/cmake.z_parT
  cmake_files/cmake.z_shift
  cmake_files/cmake.z_pseudo_deriv
  cmake_files/cmake.z_poisson
  cmake_files/cmake.z_cut_sub
  cmake_files/cmake.z_fixed_point_matt
  cmake_files/cmake.z_track_map_code
  cmake_files/cmake.z_track_da_tpsa
  cmake_files/cmake.z_matt
  cmake_files/cmake.z_simple_operations
  cmake_files/cmake.normal_simple
  cmake_files/cmake.z_vf_Lie_sub
  cmake_files/cmake.z_coast_matt_fake_maps
  cmake_files/cmake.z_radiation_matt_fake_maps
  cmake_files/cmake.z_canonize_matt_fake_maps
  cmake_files/cmake.z_resonance
  cmake_files/cmake.z_spin1
  cmake_files/cmake.z_spin_res
  cmake_files/cmake.z_benggston
)


include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
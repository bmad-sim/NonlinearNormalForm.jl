# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.8/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.8/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code

# Include any dependencies generated for this target.
include CMakeFiles/z_mono-exe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/z_mono-exe.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/z_mono-exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/z_mono-exe.dir/flags.make

CMakeFiles/z_mono-exe.dir/z_mono.f90.o: CMakeFiles/z_mono-exe.dir/flags.make
CMakeFiles/z_mono-exe.dir/z_mono.f90.o: z_mono.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/z_mono-exe.dir/z_mono.f90.o"
	/opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -Df2cFortran -DCESR_UNIX -DCESR_LINUX -u -traceback -cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none -fopenmp -DCESR_PLPLOT  -fPIC -O2  -c /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code/z_mono.f90 -o CMakeFiles/z_mono-exe.dir/z_mono.f90.o

CMakeFiles/z_mono-exe.dir/z_mono.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/z_mono-exe.dir/z_mono.f90.i"
	/opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -Df2cFortran -DCESR_UNIX -DCESR_LINUX -u -traceback -cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none -fopenmp -DCESR_PLPLOT  -fPIC -O2  -E /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code/z_mono.f90 > CMakeFiles/z_mono-exe.dir/z_mono.f90.i

CMakeFiles/z_mono-exe.dir/z_mono.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/z_mono-exe.dir/z_mono.f90.s"
	/opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -Df2cFortran -DCESR_UNIX -DCESR_LINUX -u -traceback -cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none -fopenmp -DCESR_PLPLOT  -fPIC -O2  -S /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code/z_mono.f90 -o CMakeFiles/z_mono-exe.dir/z_mono.f90.s

# Object files for target z_mono-exe
z_mono__exe_OBJECTS = \
"CMakeFiles/z_mono-exe.dir/z_mono.f90.o"

# External object files for target z_mono-exe
z_mono__exe_EXTERNAL_OBJECTS =

/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: CMakeFiles/z_mono-exe.dir/z_mono.f90.o
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: CMakeFiles/z_mono-exe.dir/build.make
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /Users/mgs255/Software/bmad/bmad-ecosystem/production/lib/libbsim.a
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /Users/mgs255/Software/bmad/bmad-ecosystem/production/lib/libxrlf03.a
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /Users/mgs255/Software/bmad/bmad-ecosystem/production/lib/libxrl.a
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /opt/homebrew/lib/libX11.dylib
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /opt/homebrew/lib/libXext.dylib
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /opt/homebrew/lib/libX11.dylib
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: /opt/homebrew/lib/libXext.dylib
/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono: CMakeFiles/z_mono-exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/z_mono-exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/z_mono-exe.dir/build: /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/production/bin/z_mono
.PHONY : CMakeFiles/z_mono-exe.dir/build

CMakeFiles/z_mono-exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/z_mono-exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/z_mono-exe.dir/clean

CMakeFiles/z_mono-exe.dir/depend:
	cd /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code /Users/mgs255/.julia/dev/NonlinearNormalForm/fpp-ptc-sandbox/code/CMakeFiles/z_mono-exe.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/z_mono-exe.dir/depend


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
CMAKE_COMMAND = /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/cmake-3.27.9-frqkedjliu3opxtrrowp6tmltgigaszp/bin/cmake

# The command to remove a file.
RM = /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/cmake-3.27.9-frqkedjliu3opxtrrowp6tmltgigaszp/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratch/ethanmar/tauStudy/MuColl-TauStudy/build

# Include any dependencies generated for this target.
include CMakeFiles/DefaultTauFinder.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/DefaultTauFinder.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/DefaultTauFinder.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DefaultTauFinder.dir/flags.make

CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o: CMakeFiles/DefaultTauFinder.dir/flags.make
CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o: /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/TauFinder.cc
CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o: CMakeFiles/DefaultTauFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/scratch/ethanmar/tauStudy/MuColl-TauStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o -MF CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o.d -o CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o -c /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/TauFinder.cc

CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/TauFinder.cc > CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.i

CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/TauFinder.cc -o CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.s

CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o: CMakeFiles/DefaultTauFinder.dir/flags.make
CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o: /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/EvaluateTauFinder.cc
CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o: CMakeFiles/DefaultTauFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/scratch/ethanmar/tauStudy/MuColl-TauStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o -MF CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o.d -o CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o -c /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/EvaluateTauFinder.cc

CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/EvaluateTauFinder.cc > CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.i

CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/EvaluateTauFinder.cc -o CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.s

CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o: CMakeFiles/DefaultTauFinder.dir/flags.make
CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o: /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/HelixClass.cc
CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o: CMakeFiles/DefaultTauFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/scratch/ethanmar/tauStudy/MuColl-TauStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o -MF CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o.d -o CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o -c /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/HelixClass.cc

CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/HelixClass.cc > CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.i

CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/HelixClass.cc -o CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.s

CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o: CMakeFiles/DefaultTauFinder.dir/flags.make
CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o: /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/SimpleLine.cc
CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o: CMakeFiles/DefaultTauFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/scratch/ethanmar/tauStudy/MuColl-TauStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o -MF CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o.d -o CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o -c /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/SimpleLine.cc

CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/SimpleLine.cc > CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.i

CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder/src/SimpleLine.cc -o CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.s

# Object files for target DefaultTauFinder
DefaultTauFinder_OBJECTS = \
"CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o" \
"CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o" \
"CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o" \
"CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o"

# External object files for target DefaultTauFinder
DefaultTauFinder_EXTERNAL_OBJECTS =

libDefaultTauFinder.so: CMakeFiles/DefaultTauFinder.dir/src/TauFinder.cc.o
libDefaultTauFinder.so: CMakeFiles/DefaultTauFinder.dir/src/EvaluateTauFinder.cc.o
libDefaultTauFinder.so: CMakeFiles/DefaultTauFinder.dir/src/HelixClass.cc.o
libDefaultTauFinder.so: CMakeFiles/DefaultTauFinder.dir/src/SimpleLine.cc.o
libDefaultTauFinder.so: CMakeFiles/DefaultTauFinder.dir/build.make
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/raida-1.11-zlntdvh7ywd5tqfy2uxwsdidewx6rwc6/lib/libRAIDA.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libHist.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libCore.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libImt.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libRIO.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libNet.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libHist.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libGraf.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libGraf3d.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libGpad.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libROOTDataFrame.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libTree.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libTreePlayer.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libRint.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libPostscript.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libMatrix.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libPhysics.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libMathCore.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libThread.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libMultiProc.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libROOTVecOps.so
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libMatrix.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libMathCore.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libImt.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libMultiProc.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libNet.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libRIO.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libThread.so.6.30.06
libDefaultTauFinder.so: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.3.1/root-6.30.06-ufi3jor5elnk5w677vdhr4r7bpj4evkv/lib/root/libCore.so.6.30.06
libDefaultTauFinder.so: CMakeFiles/DefaultTauFinder.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/scratch/ethanmar/tauStudy/MuColl-TauStudy/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX shared library libDefaultTauFinder.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DefaultTauFinder.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DefaultTauFinder.dir/build: libDefaultTauFinder.so
.PHONY : CMakeFiles/DefaultTauFinder.dir/build

CMakeFiles/DefaultTauFinder.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DefaultTauFinder.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DefaultTauFinder.dir/clean

CMakeFiles/DefaultTauFinder.dir/depend:
	cd /scratch/ethanmar/tauStudy/MuColl-TauStudy/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder /scratch/ethanmar/tauStudy/MuColl-TauStudy/DefaultTauFinder /scratch/ethanmar/tauStudy/MuColl-TauStudy/build /scratch/ethanmar/tauStudy/MuColl-TauStudy/build /scratch/ethanmar/tauStudy/MuColl-TauStudy/build/CMakeFiles/DefaultTauFinder.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/DefaultTauFinder.dir/depend


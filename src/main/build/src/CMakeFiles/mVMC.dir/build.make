# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/yoshimi/program/mVMC/src/main

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/yoshimi/program/mVMC/src/main/build

# Include any dependencies generated for this target.
include src/CMakeFiles/mVMC.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/mVMC.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/mVMC.dir/flags.make

src/CMakeFiles/mVMC.dir/vmcmain.c.o: src/CMakeFiles/mVMC.dir/flags.make
src/CMakeFiles/mVMC.dir/vmcmain.c.o: ../src/vmcmain.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yoshimi/program/mVMC/src/main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/mVMC.dir/vmcmain.c.o"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && /opt/intel/composer_xe_2015.3.187/bin/intel64/icc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/mVMC.dir/vmcmain.c.o   -c /Users/yoshimi/program/mVMC/src/main/src/vmcmain.c

src/CMakeFiles/mVMC.dir/vmcmain.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mVMC.dir/vmcmain.c.i"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && /opt/intel/composer_xe_2015.3.187/bin/intel64/icc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/yoshimi/program/mVMC/src/main/src/vmcmain.c > CMakeFiles/mVMC.dir/vmcmain.c.i

src/CMakeFiles/mVMC.dir/vmcmain.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mVMC.dir/vmcmain.c.s"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && /opt/intel/composer_xe_2015.3.187/bin/intel64/icc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/yoshimi/program/mVMC/src/main/src/vmcmain.c -o CMakeFiles/mVMC.dir/vmcmain.c.s

src/CMakeFiles/mVMC.dir/vmcmain.c.o.requires:

.PHONY : src/CMakeFiles/mVMC.dir/vmcmain.c.o.requires

src/CMakeFiles/mVMC.dir/vmcmain.c.o.provides: src/CMakeFiles/mVMC.dir/vmcmain.c.o.requires
	$(MAKE) -f src/CMakeFiles/mVMC.dir/build.make src/CMakeFiles/mVMC.dir/vmcmain.c.o.provides.build
.PHONY : src/CMakeFiles/mVMC.dir/vmcmain.c.o.provides

src/CMakeFiles/mVMC.dir/vmcmain.c.o.provides.build: src/CMakeFiles/mVMC.dir/vmcmain.c.o


src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o: src/CMakeFiles/mVMC.dir/flags.make
src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o: ../src/sfmt/SFMT.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yoshimi/program/mVMC/src/main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && /opt/intel/composer_xe_2015.3.187/bin/intel64/icc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/mVMC.dir/sfmt/SFMT.c.o   -c /Users/yoshimi/program/mVMC/src/main/src/sfmt/SFMT.c

src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mVMC.dir/sfmt/SFMT.c.i"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && /opt/intel/composer_xe_2015.3.187/bin/intel64/icc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/yoshimi/program/mVMC/src/main/src/sfmt/SFMT.c > CMakeFiles/mVMC.dir/sfmt/SFMT.c.i

src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mVMC.dir/sfmt/SFMT.c.s"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && /opt/intel/composer_xe_2015.3.187/bin/intel64/icc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/yoshimi/program/mVMC/src/main/src/sfmt/SFMT.c -o CMakeFiles/mVMC.dir/sfmt/SFMT.c.s

src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.requires:

.PHONY : src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.requires

src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.provides: src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.requires
	$(MAKE) -f src/CMakeFiles/mVMC.dir/build.make src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.provides.build
.PHONY : src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.provides

src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.provides.build: src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o


# Object files for target mVMC
mVMC_OBJECTS = \
"CMakeFiles/mVMC.dir/vmcmain.c.o" \
"CMakeFiles/mVMC.dir/sfmt/SFMT.c.o"

# External object files for target mVMC
mVMC_EXTERNAL_OBJECTS =

src/mVMC: src/CMakeFiles/mVMC.dir/vmcmain.c.o
src/mVMC: src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o
src/mVMC: src/CMakeFiles/mVMC.dir/build.make
src/mVMC: /Users/yoshimi/anaconda/lib/libmkl_intel_lp64.dylib
src/mVMC: /Users/yoshimi/anaconda/lib/libmkl_intel_thread.dylib
src/mVMC: /Users/yoshimi/anaconda/lib/libmkl_core.dylib
src/mVMC: /Users/yoshimi/anaconda/lib/libiomp5.dylib
src/mVMC: /opt/local/lib/openmpi-gcc49/libmpi.dylib
src/mVMC: src/CMakeFiles/mVMC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/yoshimi/program/mVMC/src/main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable mVMC"
	cd /Users/yoshimi/program/mVMC/src/main/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mVMC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/mVMC.dir/build: src/mVMC

.PHONY : src/CMakeFiles/mVMC.dir/build

src/CMakeFiles/mVMC.dir/requires: src/CMakeFiles/mVMC.dir/vmcmain.c.o.requires
src/CMakeFiles/mVMC.dir/requires: src/CMakeFiles/mVMC.dir/sfmt/SFMT.c.o.requires

.PHONY : src/CMakeFiles/mVMC.dir/requires

src/CMakeFiles/mVMC.dir/clean:
	cd /Users/yoshimi/program/mVMC/src/main/build/src && $(CMAKE_COMMAND) -P CMakeFiles/mVMC.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/mVMC.dir/clean

src/CMakeFiles/mVMC.dir/depend:
	cd /Users/yoshimi/program/mVMC/src/main/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/yoshimi/program/mVMC/src/main /Users/yoshimi/program/mVMC/src/main/src /Users/yoshimi/program/mVMC/src/main/build /Users/yoshimi/program/mVMC/src/main/build/src /Users/yoshimi/program/mVMC/src/main/build/src/CMakeFiles/mVMC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/mVMC.dir/depend

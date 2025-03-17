# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/f/Projects/x264/VCA-linux

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/f/Projects/x264/VCA-linux

# Include any dependencies generated for this target.
include source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/depend.make

# Include the progress variables for this target.
include source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/progress.make

# Include the compile flags for this target's objects.
include source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/flags.make

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct8.asm.o: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/flags.make
source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct8.asm.o: source/lib/analyzer/simd/dct8.asm
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/f/Projects/x264/VCA-linux/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building ASM_NASM object source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct8.asm.o"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && /usr/bin/nasm $(ASM_NASM_INCLUDES) $(ASM_NASM_FLAGS) -f elf64 -o CMakeFiles/vcaLibSimd12bit.dir/dct8.asm.o /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/dct8.asm

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/const-a.asm.o: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/flags.make
source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/const-a.asm.o: source/lib/analyzer/simd/const-a.asm
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/f/Projects/x264/VCA-linux/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building ASM_NASM object source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/const-a.asm.o"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && /usr/bin/nasm $(ASM_NASM_INCLUDES) $(ASM_NASM_FLAGS) -f elf64 -o CMakeFiles/vcaLibSimd12bit.dir/const-a.asm.o /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/const-a.asm

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/cpu-a.asm.o: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/flags.make
source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/cpu-a.asm.o: source/lib/analyzer/simd/cpu-a.asm
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/f/Projects/x264/VCA-linux/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building ASM_NASM object source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/cpu-a.asm.o"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && /usr/bin/nasm $(ASM_NASM_INCLUDES) $(ASM_NASM_FLAGS) -f elf64 -o CMakeFiles/vcaLibSimd12bit.dir/cpu-a.asm.o /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/cpu-a.asm

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.o: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/flags.make
source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.o: source/lib/analyzer/simd/dct-ssse3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/f/Projects/x264/VCA-linux/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.o"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -mssse3 -o CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.o -c /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/dct-ssse3.cpp

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.i"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -mssse3 -E /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/dct-ssse3.cpp > CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.i

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.s"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -mssse3 -S /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/dct-ssse3.cpp -o CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.s

# Object files for target vcaLibSimd12bit
vcaLibSimd12bit_OBJECTS = \
"CMakeFiles/vcaLibSimd12bit.dir/dct8.asm.o" \
"CMakeFiles/vcaLibSimd12bit.dir/const-a.asm.o" \
"CMakeFiles/vcaLibSimd12bit.dir/cpu-a.asm.o" \
"CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.o"

# External object files for target vcaLibSimd12bit
vcaLibSimd12bit_EXTERNAL_OBJECTS =

source/lib/analyzer/simd/libvcaLibSimd12bit.a: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct8.asm.o
source/lib/analyzer/simd/libvcaLibSimd12bit.a: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/const-a.asm.o
source/lib/analyzer/simd/libvcaLibSimd12bit.a: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/cpu-a.asm.o
source/lib/analyzer/simd/libvcaLibSimd12bit.a: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/dct-ssse3.cpp.o
source/lib/analyzer/simd/libvcaLibSimd12bit.a: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/build.make
source/lib/analyzer/simd/libvcaLibSimd12bit.a: source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/f/Projects/x264/VCA-linux/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libvcaLibSimd12bit.a"
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && $(CMAKE_COMMAND) -P CMakeFiles/vcaLibSimd12bit.dir/cmake_clean_target.cmake
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vcaLibSimd12bit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/build: source/lib/analyzer/simd/libvcaLibSimd12bit.a

.PHONY : source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/build

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/clean:
	cd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd && $(CMAKE_COMMAND) -P CMakeFiles/vcaLibSimd12bit.dir/cmake_clean.cmake
.PHONY : source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/clean

source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/depend:
	cd /mnt/f/Projects/x264/VCA-linux && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/f/Projects/x264/VCA-linux /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd /mnt/f/Projects/x264/VCA-linux /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd /mnt/f/Projects/x264/VCA-linux/source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/lib/analyzer/simd/CMakeFiles/vcaLibSimd12bit.dir/depend


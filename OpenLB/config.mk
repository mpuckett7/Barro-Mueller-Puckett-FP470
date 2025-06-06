# OpenLB build configuration
#
# THIS IS THE OPENMP CONFIG FILE BUILT FOR CS470 Final Project
# TO USE: Copy contents to the file named "config.mk" in
# ../openLB_original and ../OpenLB

# DO NOT CHANGE THE NAME OF config.mk, only replace the contents

# Compiler to use for C++ files, change to `mpic++` when using OpenMPI and GCC
CXX             := mpic++
# Compiler to use for C files (used for emebedded dependencies)
CC              := mpicc

# Suggested optimized build flags for GCC, consult `config/` for further examples
CXXFLAGS        := -O3 -Wall -march=native -mtune=native
# Uncomment to add debug symbols and enable runtime asserts
#CXXFLAGS        += -g -DOLB_DEBUG

# OpenLB requires support for C++17
# works in:
#  * gcc 9 or later      (https://gcc.gnu.org/projects/cxx-status.html#cxx17)
#  * icc 19.0 or later   (https://software.intel.com/en-us/articles/c17-features-supported-by-intel-c-compiler)
#  * clang 7 or later  (https://clang.llvm.org/cxx_status.html#cxx17)
CXXFLAGS        += -std=c++17

# optional linker flags
LDFLAGS         :=

# Parallelization mode, must be one of: OFF, MPI, OMP, HYBRID
# Note that for MPI and HYBRID the compiler also needs to be adapted.
# See e.g. `config/cpu_gcc_openmpi.mk`
PARALLEL_MODE   := HYBRID

# optional MPI and OpenMP flags
MPIFLAGS        := 
OMPFLAGS        := -fopenmp

# Options: CPU_SISD, CPU_SIMD, GPU_CUDA
# Both CPU_SIMD and GPU_CUDA require system-specific adjustment of compiler flags.
# See e.g. `config/cpu_simd_intel_mpi.mk` or `config/gpu_only.mk` for examples.
# CPU_SISD must always be present.
PLATFORMS       := CPU_SISD CPU_SIMD

# Fundamental arithmetic data type
# Common options are float or double
FLOATING_POINT_TYPE := double

# Any entries are passed to the compiler as `-DFEATURE_*` declarations
# Used to enable some alternative code paths and dependencies
FEATURES        :=

# Set to OFF if libz and tinyxml are provided by the system (optional)
USE_EMBEDDED_DEPENDENCIES := ON
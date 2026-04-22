#
# This file is licensed under the BSD 3-Clause License.
##/*=========================================================================
#
#
# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)

# First try to find cuFFTW (CUDA FFT Wrapper) in NVIDIA HPC SDK installation
find_library (FFTW_LIBRARIES NAMES cufftw 
		HINTS $ENV{NVHPC_ROOT}/math_libs/lib64
		      $ENV{NCCL_DIR}/math_libs/lib
    NO_DEFAULT_PATH)

# If cuFFTW not found in HPC SDK, try system paths
if (NOT FFTW_LIBRARIES)
  find_library (FFTW_LIBRARIES NAMES cufftw)
endif (NOT FFTW_LIBRARIES)

# If cuFFTW still not found, fall back to regular FFTW
if (NOT FFTW_LIBRARIES)
  find_library (FFTW_LIBRARIES NAMES fftw3)
endif (NOT FFTW_LIBRARIES)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
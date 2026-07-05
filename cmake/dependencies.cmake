# dependencies.cmake — resolve PSOPT's dependencies for a normal (non-superbuild)
# build. Each dep is found first; Eigen (header-only) and ADOL-C have a
# fetch/build fallback. Heavy deps (IPOPT, CasADi) must be found here — build
# them from source with -DPSOPT_SUPERBUILD=ON, or point at a prefix
# (e.g. -DIPOPT_ROOT=/opt/claude/ipopt, -DCMAKE_PREFIX_PATH=...).
include(FetchContent)

# ---- Eigen3 (header-only; fetch if missing) ---------------------------------
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found — fetching Eigen 3.4.0")
  FetchContent_Declare(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0 GIT_SHALLOW TRUE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(Eigen3)
endif()

# ---- BLAS / LAPACK (used by IPOPT/MUMPS/HSL/CasADi/Sacado; carried by IPOPT
#      transitively, but expose the imported targets when available) ----------
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
if(BLAS_FOUND)
  message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
endif()
if(LAPACK_FOUND)
  message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
endif()

# ---- OpenMP (opt-in): threaded BLAS/HSL(ma86/ma97)/ColPack --------------------
if(PSOPT_WITH_OPENMP)
  # Apple clang has no built-in OpenMP; point FindOpenMP at MacPorts/Homebrew libomp.
  if(APPLE)
    # MacPorts installs libomp under lib/libomp and include/libomp subdirs.
    find_path(_omp_inc omp.h HINTS ${OpenMP_ROOT} /opt/local/include/libomp /opt/local/include /usr/local/include /opt/homebrew/include)
    find_library(_omp_lib NAMES omp libomp HINTS ${OpenMP_ROOT} /opt/local/lib/libomp /opt/local/lib /usr/local/lib /opt/homebrew/lib)
    if(_omp_inc AND _omp_lib)
      foreach(_l C CXX)
        set(OpenMP_${_l}_FLAGS "-Xpreprocessor -fopenmp -I${_omp_inc}" CACHE STRING "" FORCE)
        set(OpenMP_${_l}_LIB_NAMES "omp" CACHE STRING "" FORCE)
      endforeach()
      set(OpenMP_omp_LIBRARY "${_omp_lib}" CACHE FILEPATH "" FORCE)
    endif()
  endif()
  find_package(OpenMP REQUIRED)
  message(STATUS "PSOPT_WITH_OPENMP=ON — OpenMP ${OpenMP_CXX_VERSION} (${OpenMP_omp_LIBRARY})")
endif()

# ---- MPI (opt-in): distributed-memory MUMPS ---------------------------------
if(PSOPT_WITH_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "PSOPT_WITH_MPI=ON — MPI ${MPI_CXX_VERSION}")
  message(WARNING "PSOPT_WITH_MPI: MPI-coupled IPOPT/MUMPS is known to crash at "
                  "load on macOS 26 (dyld TLS-in-constructor). Use only where verified.")
endif()

# ---- ADOL-C (custom Findadolc has its own download+build fallback) ----------
find_package(adolc REQUIRED)

# ---- IPOPT (required, PUBLIC) -----------------------------------------------
find_package(IPOPT)
if(NOT IPOPT_FOUND)
  message(FATAL_ERROR
    "IPOPT not found. Either:\n"
    "  * point at an install:  -DIPOPT_ROOT=/opt/claude/ipopt  (or add its lib/pkgconfig to PKG_CONFIG_PATH), or\n"
    "  * build everything from source:  cmake -DPSOPT_SUPERBUILD=ON ...")
endif()
if(IPOPT_VERSION)
  message(STATUS "Found IPOPT ${IPOPT_VERSION}")
endif()

# ---- HSL: lives inside IPOPT; just steer the default solver -----------------
if(PSOPT_WITH_HSL)
  if(PSOPT_DEFAULT_LINEAR_SOLVER STREQUAL "mumps")
    set(PSOPT_DEFAULT_LINEAR_SOLVER "ma97" CACHE STRING "" FORCE)
    message(STATUS "PSOPT_WITH_HSL=ON — default linear solver set to ma97 (OpenMP-parallel)")
  endif()
  message(STATUS "PSOPT_WITH_HSL=ON — the linked IPOPT must have been built with HSL "
                 "(use -DPSOPT_SUPERBUILD=ON -DCOINHSL_SOURCE_DIR=... to build it)")
endif()

# ---- CasADi (only when the CasADi backend is enabled) -----------------------
if(PSOPT_WITH_CASADI)
  find_package(casadi QUIET)
  if(NOT casadi_FOUND)
    message(FATAL_ERROR
      "PSOPT_WITH_CASADI=ON but CasADi was not found. Build it from source with "
      "-DPSOPT_SUPERBUILD=ON, or install CasADi and set casadi_DIR.")
  endif()
  message(STATUS "Found CasADi — CasADi NLP backend enabled")
endif()

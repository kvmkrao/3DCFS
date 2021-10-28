# =====================================================================
# \file CMakeLists.txt
# \author V Kotteda 
# \date 2021-02-04
# =====================================================================

cmake_minimum_required(VERSION 3.0) 
# set paths
SET(MPI_DIR /usr/) 

# set default c++ compilers
SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)
SET(CMAKE_C_COMPILER ${MPI_DIR}/bin/mpicc)

# set default flags for release and debug mode
# specify the C++ standard
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED True)
enable_language (Fortran)
# =====================================================================
set(CMAKE_Fortran_FLAGS "-g -O3") 
# Project settings
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Begin configuration")
MESSAGE("====================================")

# start generator
project(3dcfs)

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

# other cmake setting
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

#find_package(MPI REQUIRED)
# =====================================================================
#SET src directory path 
SET(SRC ${CMAKE_SOURCE_DIR}/src)

#set list of files to compile 
set(SOURCE_EXE
     ${SRC}/archinfo.F
     ${SRC}/baldwin.F
     ${SRC}/blas.F
     ${SRC}/blk1_iter_sd.F
     ${SRC}/blk1_iter_sd_mf.F
     ${SRC}/blkmvp_turb_sd.F
     ${SRC}/blkrhs_turb_sd.F
     ${SRC}/blk_turb.F
     ${SRC}/ckvol.F
     ${SRC}/csr.c
     ${SRC}/diskout_sd.F
     ${SRC}/diskout_turb_sd.F
     ${SRC}/echoinput.F
     ${SRC}/error.F
     ${SRC}/f2snew.F
     ${SRC}/fgmres_smag.F
     ${SRC}/fgmres_turb_sd.F
     ${SRC}/filter.F
     ${SRC}/frifwd.F
     ${SRC}/gen.F
     ${SRC}/gen_turb.F
     ${SRC}/gquad.F
     ${SRC}/hypo2k.F
     ${SRC}/initialize.F
     ${SRC}/keywr.F
     ${SRC}/main.F
     ${SRC}/mkshlocal.F
     ${SRC}/mylen.F
     ${SRC}/mymalloc.c
     ${SRC}/Ncomm.F
     ${SRC}/norm.F
     ${SRC}/parseinput.F
     ${SRC}/pn.F
     ${SRC}/probe.F
     ${SRC}/resout.F
     ${SRC}/rotation.F
     ${SRC}/shape1.F
     ${SRC}/shape2.F
     ${SRC}/shape3.F
     ${SRC}/shape.F
     ${SRC}/test.F
     ${SRC}/trans_data.F
     ${SRC}/transform.F
     ${SRC}/transition.F
     ${SRC}/tree.F
     ${SRC}/trip.F
     ${SRC}/tstamp.F
     ${SRC}/updated.F
     ${SRC}/utility1.c
     ${SRC}/utility.c
     ${SRC}/yon.F ) 

include_directories(".")
include_directories(/usr/local/include
                   ${MPI_DIR}/include
                   src) 

set(PROJECT_LINK_LIBS libmpi.so)
LINK_DIRECTORIES(${MPI_DIR}/lib/)

ADD_EXECUTABLE(3dcfs ${SOURCE_EXE})

TARGET_LINK_LIBRARIES(3dcfs
            -lmpi)

# =====================================================================
# Print information
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Config Information:")
MESSAGE("====================================")
MESSAGE("")

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})
MESSAGE("Build shared library: " ${BUILD_SHARED_LIBS})
MESSAGE("Doxygen available: " ${DOXYGEN_FOUND})
MESSAGE("")

MESSAGE("Compilers:")
MESSAGE("\tFortran compiler: " ${CMAKE_Fortran_COMPILER})
MESSAGE("")



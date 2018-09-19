################################################################################
# Required cmake version
################################################################################

cmake_minimum_required(VERSION 3.3.0)


################################################################################
# Project
################################################################################

project(pdp3)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# # set helper pathes to find libraries and packages
# # Add specific hints
# list(APPEND CMAKE_PREFIX_PATH "$ENV{MPI_ROOT}")
# list(APPEND CMAKE_PREFIX_PATH "$ENV{CUDA_ROOT}")
# list(APPEND CMAKE_PREFIX_PATH "$ENV{BOOST_ROOT}")
# list(APPEND CMAKE_PREFIX_PATH "$ENV{HDF5_ROOT}")
# list(APPEND CMAKE_PREFIX_PATH "$ENV{ADIOS_ROOT}")
# # Add from environment after specific env vars
# list(APPEND CMAKE_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}")


###############################################################################
# Language Flags
###############################################################################

# enforce C++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)


################################################################################
# Add include directories
################################################################################

add_subdirectory(lib/tinyxml2)
add_subdirectory(lib/tinyvec3d)

include_directories(lib/tinyxml2 lib/tinyvec3d include)


################################################################################
# Find HDF5
################################################################################

find_package(HDF5 REQUIRED COMPONENTS C)
if(HDF5_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HDF5_CXX_FLAGS}")
endif()

################################################################################
# Find OpenMP
################################################################################

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



# ################################################################################
# # Find tinyxml2
# ################################################################################

# find_package(tinyxml2)
# if(TINYXML2_FOUND)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TINYXML2_CXX_FLAGS}")
# endif()


################################################################################
# Find math from stdlib
################################################################################

if(NOT WIN32)
    # automatically added on windows
    set(LIBS ${LIBS} m)
endif()


################################################################################
# PIConGPU options
################################################################################

set(PIC_VERBOSE "1" CACHE STRING
    "Set verbosity level for PIConGPU (default is only physics output)")
add_definitions(-DPIC_VERBOSE_LVL=${PIC_VERBOSE})



################################################################################
# Warnings
################################################################################

# GNU
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
# ICC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
# PGI
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minform=inform")
endif()


################################################################################
# Compile & Link PIConGPU
################################################################################

# file(GLOB_RECURSE CUDASRCFILES "*.cu")
file(GLOB_RECURSE SRCFILES "src/*.cpp")

# add_library(picongpu-hostonly
#     STATIC
#     ${SRCFILES}
# )
# target_link_libraries(picongpu-hostonly ${LIBS})

# if("${PMACC_CUDA_COMPILER}" STREQUAL "clang")
#     add_executable(picongpu
#         ${CUDASRCFILES}
#     )

#     set_target_properties(picongpu PROPERTIES COMPILE_FLAGS ${CLANG_BUILD_FLAGS})
#     set_target_properties(picongpu PROPERTIES LINKER_LANGUAGE CXX)
#     set_source_files_properties(${CUDASRCFILES} PROPERTIES LANGUAGE CXX)

#     target_link_libraries(picongpu ${LIBS} picongpu-hostonly)
# else()


add_executable(pdp3
	${SRCFILES}
	)

target_link_libraries(pdp3 PRIVATE tinyxml2)
target_link_libraries(pdp3 PRIVATE tinyvec3d)

target_include_directories(pdp3 PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(pdp3 PRIVATE ${HDF5_C_LIBRARIES})

#     target_link_libraries(picongpu ${LIBS} picongpu-hostonly)
# endif()

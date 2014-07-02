cmake_minimum_required(VERSION 2.8)

project(DislScatHexMCMisfitThreading)

#Check if important libraries are available.

#Gnu Scientific library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules") # findGSL module is not shipped with cmake
find_package(GSL REQUIRED)

#BLAS
find_package(BLAS REQUIRED)

#MPI
find_package(MPI REQUIRED)

#BOOST
find_package(Boost COMPONENTS filesystem system REQUIRED)

#some libraries are build as external projects
include(ExternalProject)

#config++ library
set(config_project_prefix "${CMAKE_BINARY_DIR}")
ExternalProject_Add(
    config_project
    PREFIX ${config_project_prefix}
    URL http://www.hyperrealm.com/libconfig/libconfig-1.4.9.tar.gz
    SOURCE_DIR ${config_project_prefix}/src/config_project
    CONFIGURE_COMMAND ${config_project_prefix}/src/config_project/configure --prefix=${config_project_prefix}/src/config_project-build
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
)

#Miller indices library
ExternalProject_Add(
    miller_project
    PREFIX ${CMAKE_BINARY_DIR}
    GIT_REPOSITORY https://github.com/vifactor/MillerIndices.git
    GIT_TAG origin/master
)

#compile project itself
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src/config_project-build/include)
include_directories(${CMAKE_BINARY_DIR}/src/miller_project-build)

link_directories(${CMAKE_BINARY_DIR}/src/config_project-build/lib)
link_directories(${CMAKE_BINARY_DIR}/src/miller_project-build)

file(GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
set(EXEC_NAME DislScatHexMCMisfitThreading)
add_executable(${EXEC_NAME} ${CPP_SOURCES})

add_dependencies(${EXEC_NAME} config_project)
add_dependencies(${EXEC_NAME} miller_project)

target_link_libraries(${EXEC_NAME} config++ gsl miller
                                    ${Boost_LIBRARIES} 
                                    ${BLAS_LIBRARIES}
                                    ${MPI_LIBRARIES})

cmake_minimum_required(VERSION 2.8)

#Bring the headers, such as Student.h into the project
#include_directories(include)

project(DVR Fortran)
#project(DVR Fortran C)

#set(CMAKE_BUILD_TYPE Debug)

find_package(BLAS)
if (NOT BLAS_LIBRARIES)
  message(FATAL_ERROR "A required library with BLAS API not found.")
else()
  message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
endif()

set(${PROJECT_NAME}_CONFIGURATIONS dvr)

add_custom_target(get_date
  COMMAND ${PROJECT_SOURCE_DIR}/src/date.sh ${CMAKE_Fortran_COMPILER} ${CMAKE_Fortran_FLAGS}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src
  COMMENT "Getting information about the code"
)

set(CMAKE_C_FLAGS "-Wall ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS "-warn -O5")
#set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -check bounds -O5")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -check all -fpe0 -warn -traceback -debug extended -O5")
set(CMAKE_EXE_LINKER_FLAGS "-qopenmp -lblas -llapack")
#set(CMAKE_EXE_LINKER_FLAGS "-qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -O5")

file(GLOB SOURCES "src/*.f90")

add_executable(dvr ${SOURCES})
add_dependencies(dvr get_date)
target_link_libraries(dvr ${PROJECT_LINK_LIBS} )

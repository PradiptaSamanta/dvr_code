cmake_minimum_required(VERSION 2.8.4)

#Bring the headers, such as Student.h into the project
#include_directories(include)

#project(arpack Fortran)
project(DVR Fortran)

set(CMAKE_BUILD_TYPE Release)

set(${PROJECT_NAME}_CONFIGURATIONS dvr)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/build/mods/)
set(CMAKE_BINARY_DIR ${PROJECT_SOURCE_DIR}/build/bin/)

option(USE_OPENMP "Use OpenMP multi-threading" ON) 
#set(${PROJECT_NAME}_CONFIGURATIONS arpack)

#find_package(LAPACK REQUIRED)
#find_path(BLAS_INCLUDE_DIRS libblas.so
#/usr/include
#/usr/lib64
#/usr/local/include
#$ENV{BLAS_HOME}/include)

#set (BLAS_DIR "/usr/lib64/")
#find_package(BLAS REQUIRED)

#HAVE_OPENMP

set (CMAKE_FORTRAN_FLAGS "-qopenmp -lblas -llapack -O3")
set (CMAKE_EXE_LINKER_FLAGS "-qopenmp -lblas -llapack -O3")
#set (CMAKE_FORTRAN_FLAGS "-lblas -llapack -O3")
#set (CMAKE_EXE_LINKER_FLAGS "-lblas -llapack -O3")

#if (OPENMP_FOUND)
#set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} ${OpenMP_FORTRAN_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "arpack/*.f")
#include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(arpack STATIC ${SOURCES})

#install(TARGETS arpack DESTINATION ${PROJECT_SOURCE_DIR}/arpack/)
#add_executable(arpack ${SOURCES})

set (PROJECT_LINK_LIBS ${PROJECT_SOURCE_DIR}/build/libarpack.a -lblas -llapack)

# The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.f90")
#include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(dvr ${SOURCES})
target_link_libraries(dvr ${PROJECT_LINK_LIBS} )

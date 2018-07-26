cmake_minimum_required(VERSION 2.8.4)

#Bring the headers, such as Student.h into the project
#include_directories(include)

project(arpack Fortran)
set(CMAKE_BUILD_TYPE Release)

set(${PROJECT_NAME}_CONFIGURATIONS arpack)

#find_package(LAPACK REQUIRED)
#find_path(BLAS_INCLUDE_DIRS libblas.so
#/usr/include
#/usr/lib64
#/usr/local/include
#$ENV{BLAS_HOME}/include)

#set (BLAS_DIR "/usr/lib64/")
#find_package(BLAS REQUIRED)

set (CMAKE_FORTRAN_FLAGS "-lblas -llapack")
# The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "arpack/*.f")
#include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(arpack STATIC ${SOURCES})

install(TARGETS arpack DESTINATION ${PROJECT_SOURCE_DIR}/arpack/)
#add_executable(arpack ${SOURCES})



project(DVR Fortran)

set(${PROJECT_NAME}_CONFIGURATIONS dvr)

set (PROJECT_LINK_LIBS ${PROJECT_SOURCE_DIR}/build/libarpack.a -lblas)

# The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.f90")
#include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(dvr ${SOURCES})
target_link_libraries(dvr ${PROJECT_LINK_LIBS} )

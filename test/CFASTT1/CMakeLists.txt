cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME_STR CFASTT1)

#Bring the headers into the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
	${PROJECT_SOURCE_DIR}/include)

#The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

#Adding compiler flags
set (CMAKE_CXX_FLAGS "-DWM_DP -DNoRepository -DHyProRescource=${PROJECT_SOURCE_DIR}/resources -O0 -g3 -Wall -fpermissive -std=c++0x -fPIC -MMD -MP -Wno-unused-local-typedefs")

#Generate the shared library from the sources
add_executable(${PROJECT_NAME_STR} ${SOURCES})

#OpenFOAM include libraries
include_directories($ENV{FOAM_SRC}/OpenFOAM/lnInclude
	$ENV{FOAM_SRC}/finiteVolume/lnInclude
	$ENV{FOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude
	$ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
	$ENV{FOAM_SRC}/thermophysicalModels/specie/lnInclude)

#Link the libraries
target_link_libraries(${PROJECT_NAME_STR}
	HyPro
	"$ENV{FOAM_LIBBIN}/libspecie.so"
	"$ENV{FOAM_LIBBIN}/libthermophysicalFunctions.so"
	"$ENV{FOAM_LIBBIN}/libOpenFOAM.so"
	"$ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so")

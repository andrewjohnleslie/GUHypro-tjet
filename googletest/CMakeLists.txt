cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME_STR googletest)

#Bring the headers into the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/
        ${PROJECT_SOURCE_DIR}/test/Hyperion/
#         ${PROJECT_SOURCE_DIR}/test/CFD/
#         ${PROJECT_SOURCE_DIR}/test/SCRAMSPACE/
#         ${PROJECT_SOURCE_DIR}/test/STATALTEXRamjet/
#         ${PROJECT_SOURCE_DIR}/test/GP/
		)

#The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

#Disable warning regarding deprecated C++ classes.
# @TODO modify the gpc lybrary in order to avoid this warning
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -fPIC")

#Adding definitions
add_definitions(-DHyProRescource=${PROJECT_SOURCE_DIR}/resources -DSTXRescource=${PROJECT_SOURCE_DIR}/test/STATALTEXRamjet/resources -DGPRescource=${PROJECT_SOURCE_DIR}/test/GP/resources)

#Generate the executable from the sources
add_executable(${PROJECT_NAME_STR} ${SOURCES})

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

##OpenFOAM include libraries
#include_directories($ENV{FOAM_SRC}/OpenFOAM/lnInclude
#	$ENV{FOAM_SRC}/finiteVolume/lnInclude
#	$ENV{FOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude
#	$ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
#	$ENV{FOAM_SRC}/thermophysicalModels/specie/lnInclude)
#
#add_library(specie SHARED IMPORTED)
#set_property(TARGET specie PROPERTY IMPORTED_LOCATION $ENV{FOAM_LIBBIN}/libspecie.so)
#add_library(OpenFOAM SHARED IMPORTED)
#set_property(TARGET OpenFOAM PROPERTY IMPORTED_LOCATION $ENV{FOAM_LIBBIN}/libOpenFOAM.so)
#
#add_dependencies(OpenFOAM specie)

#Link the libraries
target_link_libraries(${PROJECT_NAME_STR}
	HyPro
	${GTEST_LIBRARIES}
	pthread
	boost_serialization)
#	specie
#	OpenFOAM
#	"$ENV{FOAM_LIBBIN}/libthermophysicalFunctions.so"
#	"$ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so")
	
#Link libraries
link_directories(${GTKMM_LIBRARY_DIRS})
target_link_libraries(${PROJECT_NAME_STR} ${GTKMM_LIBRARIES})

link_directories(${GOOCANVASMM_LIBRARY_DIRS})
target_link_libraries(${PROJECT_NAME_STR} ${GOOCANVASMM_LIBRARIES})

add_library(GPC SHARED IMPORTED)
set_property(TARGET GPC PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/ThirdParty/gpc/lib/libgp.so)
target_link_libraries(${PROJECT_NAME_STR} GPC)

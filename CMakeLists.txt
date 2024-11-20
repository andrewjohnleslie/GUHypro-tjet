cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME_STR HyPro)
project(${PROJECT_NAME_STR})
# The version number.
set (HyPro_VERSION_MAJOR 2)
set (HyPro_VERSION_MINOR 0)
message(--------------------------------------)
# Set CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#Allows optional building of matlab executables
option(matlab "Build the .mexa64 executables." OFF)
option(python "Build the python exectuables." OFF)
option(gtest "Build the google tests." ON)
message("Build MATLAB: " ${matlab})
message("Build Python: " ${python})
message("Build Google Tests: " ${gtest})

#Update symbolic links in include directory
file(GLOB INCLUDE ${PROJECT_SOURCE_DIR}/src/*.h
	${PROJECT_SOURCE_DIR}/src/inletModels/*.h
	${PROJECT_SOURCE_DIR}/src/nozzleModels/*.h
	${PROJECT_SOURCE_DIR}/src/numerics/*.h
	${PROJECT_SOURCE_DIR}/src/thermo/*.h
	${PROJECT_SOURCE_DIR}/src/optimization/*.h
	${PROJECT_SOURCE_DIR}/src/Graphics/*.h
	${PROJECT_SOURCE_DIR}/src/combustion/*.h
	${PROJECT_SOURCE_DIR}/src/core/*.h
	${PROJECT_SOURCE_DIR}/src/injectors/*.h
	${PROJECT_SOURCE_DIR}/src/solvers/*.h
	${PROJECT_SOURCE_DIR}/src/utilities/*.h
	${PROJECT_SOURCE_DIR}/src/trajectoryMechanics/*.h
    ${PROJECT_SOURCE_DIR}/src/gasTurbine/*.h)

#Bring the headers into the project
include_directories(${PROJECT_SOURCE_DIR}/src/
		${PROJECT_SOURCE_DIR}/src/inletModels
		${PROJECT_SOURCE_DIR}/src/nozzleModels
		${PROJECT_SOURCE_DIR}/src/numerics
		${PROJECT_SOURCE_DIR}/src/thermo
		${PROJECT_SOURCE_DIR}/src/optimization
		${PROJECT_SOURCE_DIR}/src/Graphics
		${PROJECT_SOURCE_DIR}/src/combustion
		${PROJECT_SOURCE_DIR}/src/core
		${PROJECT_SOURCE_DIR}/src/injectors
		${PROJECT_SOURCE_DIR}/src/solvers
		${PROJECT_SOURCE_DIR}/src/utilities
		${PROJECT_SOURCE_DIR}/src/trajectoryMechanics
        ${PROJECT_SOURCE_DIR}/src/gasTurbine)

#The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/inletModels/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/nozzleModels/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/numerics/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/thermo/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/optimization/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/Graphics/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/Graphics/*.cc"
	"${PROJECT_SOURCE_DIR}/src/combustion/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/core/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/injectors/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/solvers/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/utilities/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/trajectoryMechanics/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/gasTurbine/*.cpp")

#Adding compiler flags
# set (CMAKE_CXX_FLAGS "-O1 -Wall -std=c++0x -fPIC -MMD -MP -Wno-unused-local-typedefs -Wno-reorder")
set (CMAKE_CXX_FLAGS "-O0 -g3 -Wall -std=c++14 -fPIC -MMD -MP -Wno-unused-local-typedefs -Wno-reorder")

#Generate the shared library from the sources
add_library(${PROJECT_NAME_STR} SHARED ${SOURCES})

#Create symbolic links between header files and include dir
if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/include)
  add_custom_command(
    TARGET HyPro PRE_BUILD
    COMMAND rm -r ${PROJECT_SOURCE_DIR}/include)
endif()
add_custom_command(
  TARGET HyPro PRE_BUILD
  COMMAND mkdir -p ${PROJECT_SOURCE_DIR}/include
  COMMAND ln -sf ${INCLUDE} ${PROJECT_SOURCE_DIR}/include)

message("-- Update or create include directory")
if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/include)
  message("-- Update or create include directory -- update")
else()
  message("-- Update or create include directory -- create")
endif()

find_package (CANTERA)
if (CANTERA_FOUND)
	include_directories(${CANTERA_INCLUDE_DIRS})
	target_link_libraries(${PROJECT_NAME_STR} ${CANTERA_LIBRARIES})
endif (CANTERA_FOUND)

find_package (Boost COMPONENTS system serialization REQUIRED)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME_STR} ${Boost_LIBRARIES})
endif (Boost_FOUND)

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
	include_directories(${LAPACK_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME_STR} ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)

find_package(Sundials)
if (SUNDIALS_FOUND)
	include_directories(${SUNDIALS_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME_STR} ${SUNDIALS_LIBRARIES})
endif (SUNDIALS_FOUND)

# find_package(Threads REQUIRED)
# if(THREADS_HAVE_PTHREAD_ARG)
# 	target_compile_options(PUBLIC {PROJECT_NAME_STR} "-pthread")
# endif()
# if(CMAKE_THREAD_LIBS_INIT)
# 	target_link_libraries(${PROJECT_NAME_STR} "${CMAKE_THREAD_LIBS_INIT}")
# endif()
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME_STR} Threads::Threads)


add_definitions(-DWM_DP -DNoRepository)
add_definitions("-DHyProRescource=${PROJECT_SOURCE_DIR}/resources")

# #PLEASE NOTE -- IN ORDER TO BUILD MATLAB WRAPPERS:
# #LIBHYPRO.SO MUST BE ON PATH
# #MEX MUST BE ON PATH

find_package(PkgConfig REQUIRED) # Check pkg-config is available. It is used for gtk and goocanvas

pkg_check_modules(GTKMM gtkmm-3.0 REQUIRED)
include_directories(${GTKMM_INCLUDE_DIRS})

pkg_check_modules(GOOCANVASMM goocanvasmm-2.0 REQUIRED)
include_directories(${GOOCANVASMM_INCLUDE_DIRS})

# Include Cairo, the pkg-config was supposed to include it but it does not work
include_directories(/usr/include/cairomm-1.0/cairomm)

# Include gpc library
include_directories(${PROJECT_SOURCE_DIR}/ThirdParty/gpc/include)


# Include fmt library
include_directories(${PROJECT_SOURCE_DIR}/ThirdParty/fmt-4.1.0)


# Add boost serialization
#find_package(Boost REQUIRED COMPONENTS serialization)
#target_link_libraries(${PROJECT_NAME_STR}
#        ${Boost_LIBRARIES}
#)

#Add subdirectories which depend on HyPro being built
add_subdirectory(test)

find_package(GTest)
if(GTEST_FOUND)
  if(gtest)
        add_subdirectory(googletest)
  endif()
endif()

find_package(Matlab)
if (MATLAB_FOUND)
	include_directories(${MATLAB_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME_STR} ${MATLAB_LIBRARIES})
endif (MATLAB_FOUND)

# Custom additions for finding fmt and yaml-cpp
find_package(fmt REQUIRED)
find_package(YamlCpp REQUIRED)
include_directories(${YAMLCPP_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME_STR} fmt::fmt ${YAMLCPP_LIBRARIES})

message(--------------------------------------)
	# #PLEASE NOTE -- IN ORDER TO BUILD MATLAB WRAPPERS:
# #LIBHYPRO.SO MUST BE ON PATH
# #MEX MUST BE ON PATH

#In order to build the matlab wrappers, 'make matlab' should be run after cmake
#add_custom_target(matlab)
#add_dependencies(matlab mHyperion)
#add_dependencies(matlab mCFD)
#add_dependencies(matlab mSCRAMSPACE)
#add_dependencies(matlab mSTATALTEXRamjet)

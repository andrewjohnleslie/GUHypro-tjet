cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME_STR OpenModel)

#Bring the headers into the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)

#Generate the shared library from the sources
add_executable(${PROJECT_NAME_STR} ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME_STR}.cpp)

#Link the libraries
target_link_libraries(${PROJECT_NAME_STR}
        HyPro
        specie
        OpenFOAM
        boost_serialization
        "$ENV{FOAM_LIBBIN}/libthermophysicalFunctions.so"
        "$ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so")

if(matlab)
	add_custom_target(
			${PROJECT_NAME_STR}.mexa64 ALL
			COMMAND "HyPro=${PROJECT_SOURCE_DIR}" ${CMAKE_MAKE_PROGRAM}
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/matlab
			COMMENT "Original Makefile for ${PROJECT_NAME_STR}.mexa64"
			COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/matlab/${PROJECT_NAME_STR}.mexa64 ${CMAKE_CURRENT_BINARY_DIR}
			DEPENDS HyPro
	)
endif()

target_link_libraries(${PROJECT_NAME_STR} ${GTKMM_LIBRARIES})

target_link_libraries(${PROJECT_NAME_STR} ${GOOCANVASMM_LIBRARIES})

target_link_libraries(${PROJECT_NAME_STR} GPC)


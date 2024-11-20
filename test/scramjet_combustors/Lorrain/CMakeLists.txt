cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME_STR lorrain)

#Bring the headers into the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/
        ${PROJECT_SOURCE_DIR}/include)

#Generate the shared library from the sources
add_executable(${PROJECT_NAME_STR} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/${PROJECT_NAME_STR}.cpp)

find_package (CANTERA)
message(status "** CANTERA Include: ${CANTERA_INCLUDE_DIR}")
message(status "** CANTERA Libraries: ${CANTERA_LIBRARY_DIRS}")
message(status "** CANTERA Libraries: ${CANTERA_LIBRARIES}")
if (CANTERA_FOUND)
    include_directories(${CANTERA_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME_STR} ${CANTERA_LIBRARIES})
endif (CANTERA_FOUND)

#add_dependencies(HyPro ${CANTERA_LIBRARIES})

find_package (Boost COMPONENTS system REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME_STR} ${Boost_LIBRARIES})
endif (Boost_FOUND)

message(status "** LAPACK Include: ${LAPACK_INCLUDE_DIR}")
message(status "** LAPACK Libraries: ${LAPACK_LIBRARY_DIRS}")
message(status "** LAPACK Libraries: ${LAPACK_LIBRARIES}")
find_package(LAPACK)
if (LAPACK_FOUND)
    include_directories(${LAPACK_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME_STR} ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)

find_package(Sundials)
message(status "** Sundials Include: ${SUNDIALS_INCLUDE_DIR}")
message(status "** Sundials Libraries: ${SUNDIALS_LIBRARY_DIRS}")
message(status "** Sundials Libraries: ${SUNDIALS_LIBRARIES}")
if (SUNDIALS_FOUND)
    include_directories(${SUNDIALS_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME_STR} ${SUNDIALS_LIBRARIES})
endif (SUNDIALS_FOUND)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(PUBLIC reactor_test "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(${PROJECT_NAME_STR} "${CMAKE_THREAD_LIBS_INIT}")
endif()

#Link the libraries
target_link_libraries(${PROJECT_NAME_STR}
        HyPro
        boost_serialization)

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

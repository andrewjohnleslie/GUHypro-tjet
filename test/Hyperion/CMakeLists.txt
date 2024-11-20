cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME_STR Hyperion)

#Bring the headers into the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/
        ${PROJECT_SOURCE_DIR}/include)

#Adding compiler flags
#set (CMAKE_CXX_FLAGS "-DWM_DP -DNoRepository -DHyProRescource=${PROJECT_SOURCE_DIR}/resources -O0 -g3 -Wall -fpermissive -std=c++0x -fPIC -MMD -MP -Wno-unused-local-typedefs")

#Generate the shared library from the sources
add_executable(${PROJECT_NAME_STR} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/${PROJECT_NAME_STR}.cpp)

find_package (CANTERA)
if (CANTERA_FOUND)
    include_directories(${CANTERA_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME_STR} ${CANTERA_LIBRARIES})
    MESSAGE(here)
endif (CANTERA_FOUND)

#add_dependencies(HyPro ${CANTERA_LIBRARIES})

#find_package (Boost COMPONENTS system python-py35 REQUIRED)
#if (Boost_FOUND)
#    include_directories("${Boost_INCLUDE_DIR}" "/usr/include/python3.5m")
#    target_link_libraries(${PROJECT_NAME_STR} ${Boost_LIBRARIES})
#endif (Boost_FOUND)
#
#FIND_PACKAGE(PythonLibs 3.5 REQUIRED)

find_package (Boost COMPONENTS system REQUIRED)
if (Boost_FOUND)
    include_directories("${Boost_INCLUDE_DIR}" "/usr/include/python3.5m")
    target_link_libraries(${PROJECT_NAME_STR} ${Boost_LIBRARIES})
endif (Boost_FOUND)

find_package(LAPACK)
if (LAPACK_FOUND)
    include_directories(${LAPACK_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME_STR} ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)

#find_package(Threads REQUIRED)
#if(THREADS_HAVE_PTHREAD_ARG)
#    target_compile_options(PUBLIC test_cantera "-pthread")
#endif()
#if(CMAKE_THREAD_LIBS_INIT)
#    target_link_libraries(${PROJECT_NAME_STR} "${CMAKE_THREAD_LIBS_INIT}")
#endif()
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME_STR} Threads::Threads)

find_package(Sundials)
if (SUNDIALS_FOUND)
    include_directories(${SUNDIALS_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME_STR} ${SUNDIALS_LIBRARIES})
endif (SUNDIALS_FOUND)

find_package(Matlab)
if (MATLAB_FOUND)
    include_directories(${MATLAB_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME_STR} ${MATLAB_LIBRARIES})
endif (MATLAB_FOUND)

#Link the libraries
target_link_libraries(${PROJECT_NAME_STR}
        HyPro)

target_link_libraries(${PROJECT_NAME_STR} ${GTKMM_LIBRARIES})

target_link_libraries(${PROJECT_NAME_STR} ${GOOCANVASMM_LIBRARIES})

target_link_libraries(${PROJECT_NAME_STR} GPC)

message(${matlab})
if(matlab)
#	add_custom_target(
#		${PROJECT_NAME_STR}.mexa64 ALL
#		COMMAND "HyPro=${PROJECT_SOURCE_DIR}" ${CMAKE_MAKE_PROGRAM}
#		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/matlab
#		COMMENT "Original Makefile for ${PROJECT_NAME_STR}.mexa64"
#		COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/matlab/${PROJECT_NAME_STR}.mexa64 ${CMAKE_CURRENT_BINARY_DIR}
#		DEPENDS HyPro
#	)
    matlab_add_mex(
            NAME ${PROJECT_NAME_STR}mex
            SRC ${CMAKE_CURRENT_SOURCE_DIR}/matlab/${PROJECT_NAME_STR}.cpp
            LINK_TO ${CANTERA_LIBRARIES} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${SUNDIALS_LIBRARIES} ${MATLAB_LIBRARIES} HyPro GPC ${GTKMM_LIBRARIES} ${GOOCANVASMM_LIBRARIES} boost_serialization
    )
    #configure_file(${CMAKE_CURRENT_SOURCE_DIR}/matlab/Hyperionmex.m {CMAKE_CURRENT_BINARY_DIR}/Hyperionmex.m COPYONLY)
    #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/matlab/Hyperionmex.m
    #    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_command(TARGET ${PROJECT_NAME_STR}mex POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/matlab/Hyperionmex.m $<TARGET_FILE_DIR:${PROJECT_NAME_STR}mex>)
endif()

# if(python)
#     add_library(Hyperion_python SHARED python/Hyperion.cpp)
#     target_link_libraries(Hyperion_python ${CANTERA_LIBRARIES})
#     target_link_libraries(Hyperion_python ${Boost_LIBRARIES})
#     target_link_libraries(Hyperion_python ${LAPACK_LIBRARIES})
#     target_link_libraries(Hyperion_python ${CMAKE_THREAD_LIBS_INIT})
#     target_link_libraries(Hyperion_python ${SUNDIALS_LIBRARIES})
#     target_link_libraries(Hyperion_python ${MATLAB_LIBRARIES})
#     target_link_libraries(Hyperion_python HyPro)
#     target_link_libraries(Hyperion_python GPC)
#     target_link_libraries(Hyperion_python ${GTKMM_LIBRARIES})
#     target_link_libraries(Hyperion_python ${GOOCANVASMM_LIBRARIES})
#     target_link_libraries(Hyperion_python boost_serialization)
#     target_link_libraries(Hyperion_python boost_python)
#
#     target_link_libraries(Hyperion_python ${PYTHON_LIBRARIES})
# endif()

# message(${python})
# if(python)
#     add_library(pypro SHARED pypro.cpp)
#     target_link_libraries(pypro ${CANTERA_LIBRARIES})
#     target_link_libraries(pypro ${Boost_LIBRARIES})
#     target_link_libraries(pypro ${LAPACK_LIBRARIES})
#     target_link_libraries(pypro ${CMAKE_THREAD_LIBS_INIT})
#     target_link_libraries(pypro ${SUNDIALS_LIBRARIES})
#     target_link_libraries(pypro ${MATLAB_LIBRARIES})
#     target_link_libraries(pypro HyPro)
#     target_link_libraries(pypro GPC)
#     target_link_libraries(pypro ${GTKMM_LIBRARIES})
#     target_link_libraries(pypro ${GOOCANVASMM_LIBRARIES})
#     target_link_libraries(pypro boost_serialization)
#     target_link_libraries(pypro boost_python)

#     target_link_libraries(pypro ${PYTHON_LIBRARIES})
# endif()
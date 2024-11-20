#add_library(specie SHARED IMPORTED)
#set_property(TARGET specie PROPERTY IMPORTED_LOCATION $ENV{FOAM_LIBBIN}/libspecie.so)
#add_library(OpenFOAM SHARED IMPORTED)
#set_property(TARGET OpenFOAM PROPERTY IMPORTED_LOCATION $ENV{FOAM_LIBBIN}/libOpenFOAM.so)
#
#add_dependencies(OpenFOAM specie)

#Link libraries
link_directories(${GTKMM_LIBRARY_DIRS})

link_directories(${GOOCANVASMM_LIBRARY_DIRS})

add_library(GPC SHARED IMPORTED)
set_property(TARGET GPC PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/ThirdParty/gpc/lib/libgp.so)

#Add sub directories
 add_subdirectory(HyperionExpanded)
 add_subdirectory(CFD)
 add_subdirectory(TurboJet)
 add_subdirectory(HPLPTurboJet)
 add_subdirectory(Chiara)
 add_subdirectory(TrajOpt)
 add_subdirectory(Np3)
 add_subdirectory(pyCycleTjetComparison)
 add_subdirectory(pyCycleTjetMatching)
 add_subdirectory(Enthalpy)
 add_subdirectory(TurboFan)
 
 
 # add_subdirectory(jsontest)
#add_subdirectory(HPPMM_test)
#add_subdirectory(CFASTT1)
#add_subdirectory(MarquardtRamjet)
#add_subdirectory(OpenModel)
#add_subdirectory(OpenModelHyperion)
#add_subdirectory(PipeFlow)
#add_subdirectory(Scramjet)
 add_subdirectory(SCRAMSPACE)    #free
# add_subdirectory(STATALTEXRamjet)
 add_subdirectory(GP)            #free
#add_subdirectory(Serialization)
 add_subdirectory(Rocket_test)   #free
# add_subdirectory(reactor_test)
# add_subdirectory(scramjet_combustors/Burrows)
# add_subdirectory(scramjet_combustors/Lorrain)
# add_subdirectory(scramjet_combustors/TotalT)
 add_subdirectory(Hyperion)      #free
# add_subdirectory(memleaktest)


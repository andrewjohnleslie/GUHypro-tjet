# module for Cantera, http://cantera.github.com/docs/sphinx/html/index.html
include(FindPackageHandleStandardArgs)

IF(NOT CANTERA_DIR)
   SET(CANTERA_DIR "$ENV{CANTERA_DIR}")
ENDIF()


find_path(CANTERA_INCLUDE_DIR cantera/thermo/ThermoPhase.h
  PATHS /usr/local/include/cantera #/usr/local/lib/cantera/2.0.2 /usr/local/include/cantera /usr/include/cantera
  HINTS ${CANTERA_DIR}
  PATH_SUFFIXES cantera include
)
find_library(CANTERA_LIBRARY
  NAMES cantera
  PATHS /usr/local/lib #/usr/local/lib/cantera/2.0.2/lib /usr/local/lib /usr/lib
  HINTS ${CANTERA_DIR}

  PATH_SUFFIXES cantera lib 
)

find_package_handle_standard_args(cantera  DEFAULT_MSG
                                  CANTERA_LIBRARY CANTERA_INCLUDE_DIR)


IF(CANTERA_FOUND)

  SET(CANTERA_LIBRARIES
    ${CANTERA_LIBRARY}
  )
  SET(CANTERA_INCLUDE_DIRS
    ${CANTERA_INCLUDE_DIR}
  )
  mark_as_advanced(CANTERA_INCLUDE_DIRS CANTERA_LIBRARIES )
ELSE()
  SET(CANTERA_DIR "" CACHE PATH
    "An optional hint to the cantera installation directory"
    )
  message(FATAL_ERROR "Cannot find CANTERA!")
ENDIF()

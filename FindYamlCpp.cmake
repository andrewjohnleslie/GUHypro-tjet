# Find the YAML-CPP library
find_path(YAMLCPP_INCLUDE_DIR yaml-cpp/yaml.h)
find_library(YAMLCPP_LIBRARY yaml-cpp)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YamlCpp DEFAULT_MSG YAMLCPP_LIBRARY YAMLCPP_INCLUDE_DIR)

if(YAMLCPP_FOUND)
  set(YAMLCPP_LIBRARIES ${YAMLCPP_LIBRARY})
  set(YAMLCPP_INCLUDE_DIRS ${YAMLCPP_INCLUDE_DIR})
endif()


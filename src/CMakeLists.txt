add_library(frictionless frictionless.cpp)
target_include_directories(frictionless PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

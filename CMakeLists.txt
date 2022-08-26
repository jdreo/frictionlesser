
######################################################################################
# Project settings
######################################################################################

cmake_minimum_required(VERSION 3.9)

# Set a name and a version number for your project:
project(Frictionless VERSION 0.0.1 LANGUAGES CXX)

# Initialize some default paths
include(GNUInstallDirs)

# Define the minimum C++ standard that is required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


######################################################################################
# Configurable user settings
######################################################################################

# Compilation options
option(BUILD_DOCS   "Enable building of documentation"   ON)


######################################################################################
# Dependencies
######################################################################################

# ParadisEO
set(PARADISEO_ROOT "PARADISEO_NOT_FOUND" CACHE PATH "Where to find ParadisEO")
set(PARADISEO_BUILD "${PARADISEO_ROOT}/build" CACHE PATH "Build dir of ParadisEO")

include_directories(${PARADISEO_ROOT})
include_directories(${PARADISEO_ROOT}/eo/src)
include_directories(${PARADISEO_ROOT}/mo/src)
link_directories(${PARADISEO_BUILD}/lib)
set(PARADISEO_LIBRARIES ga eoutils eo)

# Clutchlog
include_directories(ext/clutchlog)

######################################################################################
# Start building
######################################################################################

# compile the library
add_subdirectory(src)

# compile the application
add_subdirectory(app)

# compile the tests
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(ext/Catch2)
  include(./ext/Catch2/extras/Catch.cmake)
  add_subdirectory(tests)
endif()

if(BUILD_DOCS)
  # Add the documentation
  add_subdirectory(doc)
endif()

# Add an alias target for use if this project is included as a subproject in another project
add_library(frictionless::frictionless ALIAS frictionless)

# Install targets and configuration
install(
  TARGETS frictionless
  EXPORT frictionless-config
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
  EXPORT frictionless-config
  NAMESPACE frictionless::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/frictionless
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# This prints a summary of found dependencies
include(FeatureSummary)
feature_summary(WHAT ALL)


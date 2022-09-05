
######################################################################################
# Project settings
######################################################################################

cmake_minimum_required(VERSION 3.9)

# Set a name and a version number for your project:
project(Frictionless VERSION 0.0.1 LANGUAGES CXX)

# Initialize some default paths
include(GNUInstallDirs)

# Define the minimum C++ standard that is required
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


######################################################################################
# Configurable user settings
######################################################################################

# Compilation options
option(BUILD_DOCS   "Enable building of documentation"   ON)

# Common
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Weverything -Wno-c++98-compat")
# Clang
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-old-style-cast -Wno-padded -Wno-extra-semi-stmt -Wno-weak-vtables")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_CLUTCHLOG")

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

# Single-header dependencies.
include_directories(ext/clutchlog)
include_directories(ext/exceptions)

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


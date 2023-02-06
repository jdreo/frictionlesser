add_library(frictionless
    frictionless.cpp
    transcriptome.cpp
    score.cpp
    parser.cpp
    eval.cpp
    log.cpp
    cache.cpp
)
target_include_directories(frictionless PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(R
    pgamma.cpp
)
target_include_directories(R PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


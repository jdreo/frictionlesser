add_executable(tests tests.cpp frictionless_t.cpp)
target_link_libraries(tests PUBLIC frictionless Catch2::Catch2 ${PARADISEO_LIBRARIES} stdc++fs)

# allow user to run tests with `make test` or `ctest`
catch_discover_tests(tests)


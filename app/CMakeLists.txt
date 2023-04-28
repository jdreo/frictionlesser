add_executable(datatester datatester.cpp)
target_link_libraries(datatester frictionless R ${PARADISEO_LIBRARIES} stdc++fs)

add_executable(frictionlesser frictionlesser.cpp)
target_link_libraries(frictionlesser frictionless R ${PARADISEO_LIBRARIES} stdc++fs)

add_executable(frictionlesser-objfunc frictionlesser-objfunc.cpp)
target_link_libraries(frictionlesser-objfunc frictionless R ${PARADISEO_LIBRARIES} stdc++fs)

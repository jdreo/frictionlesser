add_executable(frictionlesser frictionlesser.cpp)
target_link_libraries(frictionlesser frictionless ${PARADISEO_LIBRARIES} stdc++fs)

add_executable(rank-transcriptome rank-transcriptome.cpp)
target_link_libraries(rank-transcriptome frictionless stdc++fs)

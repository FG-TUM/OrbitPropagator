set(BINARY ${PROJECT_BINARY}_test)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false
    # files in test directory
     *.h
     *.cpp
    # files in test src directory
     ../src/satellitePropagator/*.h
     ../src/satellitePropagator/*.cpp
     )

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES})

include(GoogleTest)
gtest_discover_tests(${BINARY})

find_package(heyoka)
if(heyoka_FOUND)
    target_link_libraries(${BINARY} PUBLIC ${PROJECT_BINARY}_lib gtest gmock heyoka::heyoka)
    add_compile_definitions(HEYOKA_IS_PRESENT)
else()
    target_link_libraries(${BINARY} PUBLIC ${PROJECT_BINARY}_lib gtest gmock)
endif()

find_package(Catch2 3 QUIET)
if (NOT Catch2_FOUND)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2
    GIT_TAG v3.0.0-preview3
  )

  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

set(${PROJECT_NAME}_TESTS_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cxx
)

add_executable(${PROJECT_NAME}_TESTS ${${PROJECT_NAME}_TESTS_SOURCES})
target_link_libraries(${PROJECT_NAME}_TESTS 
  PRIVATE
    ${PROJECT_NAME}
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
enable_testing()
catch_discover_tests(${PROJECT_NAME}_TESTS
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)

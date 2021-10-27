cmake_minimum_required(VERSION 3.17)
project(ashera VERSION 0.1.0 
               LANGUAGES C CXX
               DESCRIPTION "Ashera is a read correction tool")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(cmake/Warnings.cmake)
find_package(cxxopts 3.0.0 QUIET)
if (NOT cxxopts_FOUND)
  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts
    GIT_TAG v3.0.0
  )

  FetchContent_MakeAvailable(cxxopts)
endif()

find_package(bioparser 3.0.15 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.15
  ) 

  FetchContent_MakeAvailable(bioparser)
endif()

find_package(ram 2.1.1 QUIET)
if (NOT ram_FOUND)
  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/lbcb-sci/ram
    GIT_TAG 2.1.1
  )

  FetchContent_MakeAvailable(ram)
endif()

find_package(biosoup 0.10.0 QUIET)
if(NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup/releases
    GIT_TAG 0.10.0
  )

  FetchContent_MakeAvailable(biosoup)
endif()

find_package(thread_pool 4.0.0 QUIET)
if (NOT thread_pool_FOUND)
  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/rvaser/thread_pool
    GIT_TAG 4.0.0
  )
endif()

set(${PROJECT_NAME}_LIB_SRC
  src/ashera/engine.cxx
  src/ashera/io.cxx
  src/ashera/thread_pool.cxx
)

add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_LIB_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
set_default_warnings(${PROJECT_NAME} INTERFACE FALSE)

target_link_libraries(${PROJECT_NAME} 
  PUBLIC biosoup thread_pool
  PRIVATE bioparser ram)

add_executable(${PROJECT_NAME}_exe src/main.cxx)
target_link_libraries(${PROJECT_NAME}_exe PRIVATE ${PROJECT_NAME} cxxopts)

set_property(TARGET ${PROJECT_NAME}_exe PROPERTY CXX_STANDARD 17)

if (BUILD_TESTS)
  add_subdirectory(test)
endif()

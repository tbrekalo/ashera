cmake_minimum_required(VERSION 3.17)
project(ashera VERSION 0.1.0 
               LANGUAGES C CXX
               DESCRIPTION "Ashera is a read correction tool")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

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

find_package(spoa 4.0.7 QUIET)
if (NOT spoa_FOUND)
  FetchContent_Declare(
    spoa 
    GIT_REPOSITORY https://github.com/rvaser/spoa
    GIT_TAG 4.0.7
  )

  FetchContent_MakeAvailable(spoa)
endif()

find_package(biosoup 0.10.0 QUIET)
if(NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
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

  FetchContent_MakeAvailable(thread_pool)
endif()

find_package(fmt 8.0.1 QUIET)
if (NOT FMT_FOUND)
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG 8.0.1
  )

  FetchContent_MakeAvailable(fmt)
endif()

find_package(edlib 1.2.7 QUIET)
if (NOT edlib_FOUND)
  FetchContent_Declare(
    edlib
    GIT_REPOSITORY https://github.com/Martinsos/edlib
    GIT_TAG v1.2.7
  )

  FetchContent_MakeAvailable(edlib)
endif()

set(${PROJECT_NAME}_LIB_SRC
  src/ashera/detail/mapping_algo.cxx
  src/ashera/detail/overlap.cxx
  src/ashera/detail/snp_filter.cxx
  
  src/ashera/algorithm.cxx
  src/ashera/configs.cxx
  src/ashera/io.cxx
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
  PRIVATE bioparser edlib fmt spoa ram)

add_executable(${PROJECT_NAME}_exe src/main.cxx)
target_link_libraries(${PROJECT_NAME}_exe PRIVATE ${PROJECT_NAME} cxxopts fmt)

set_property(TARGET ${PROJECT_NAME}_exe PROPERTY CXX_STANDARD 17)

if (BUILD_TESTS)
  add_subdirectory(test)
endif()

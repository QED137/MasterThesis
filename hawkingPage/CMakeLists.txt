cmake_minimum_required(VERSION 3.10)
project(HawkingRadiation VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Find required libraries
find_package(Threads)

# Option for OpenMP (parallel computation)
option(USE_OPENMP "Enable OpenMP for parallel computation" OFF)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found - parallel computation enabled")
    else()
        message(WARNING "OpenMP not found - parallel computation disabled")
    endif()
endif()

# Source files
set(SOURCES
    black_hole.cpp
)

set(HEADERS
    black_hole.h
)

# Create executable
add_executable(hawking_radiation 
    ${SOURCES}
    main.cpp
)

# Link libraries
target_link_libraries(hawking_radiation PRIVATE m)

if(OpenMP_CXX_FOUND AND USE_OPENMP)
    target_link_libraries(hawking_radiation PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(hawking_radiation PRIVATE USE_OPENMP)
endif()

if(Threads_FOUND)
    target_link_libraries(hawking_radiation PRIVATE Threads::Threads)
endif()

# Set output directory
set_target_properties(hawking_radiation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Installation
install(TARGETS hawking_radiation
    RUNTIME DESTINATION bin
)

# Print configuration summary
message(STATUS "===========================================")
message(STATUS "Hawking Radiation Calculator Configuration")
message(STATUS "===========================================")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP: ${USE_OPENMP}")
message(STATUS "===========================================")
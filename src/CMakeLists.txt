cmake_minimum_required(VERSION 3.10)

# Set project name
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${ProjectId} LANGUAGES CXX)

# Debug messages
message(STATUS "-- Configuring project: ${ProjectId}")

# Gather all source files (cpp, h)
file(GLOB_RECURSE ALL_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

# Exclude CMake-related files if necessary
foreach (SOURCE_FILE ${ALL_SOURCE_FILES})
    string(FIND ${SOURCE_FILE} "/CMakeFiles/" EXCLUDE_CMAKE_FILES)
    if (NOT ${EXCLUDE_CMAKE_FILES} EQUAL -1)
        list(REMOVE_ITEM ALL_SOURCE_FILES ${SOURCE_FILE})
    endif ()
endforeach ()

# Create the executable for the non-OpenMP version
add_executable(${ProjectId}_no_omp ${ALL_SOURCE_FILES})

# If USE_OMP is set, enable OpenMP and create the OpenMP version of the executable
find_package(OpenMP REQUIRED)

# Create OpenMP executable
add_executable(${ProjectId}_omp ${ALL_SOURCE_FILES})

# Link OpenMP to the OpenMP executable
target_link_libraries(${ProjectId}_omp OpenMP::OpenMP_CXX)

# Done

cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)

project(tb2025_feb_desy)

# Set ROOT installation path explicitly
set(ROOT_DIR "/home/shihai/sw/root/root_install")  # Adjust to your ROOT installation path
list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})

# Find ROOT with correct version
find_package(ROOT REQUIRED COMPONENTS
    Core 
    RIO
    Tree
    Gpad    
    Hist
)

include_directories(${ROOT_INCLUDE_DIRS})

# Ensure FetchContent is available
include(FetchContent)

# Fetch argparse
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG        master  # Use a specific tag for more stability
)
FetchContent_MakeAvailable(argparse)

# Fetch easylogging++
FetchContent_Declare(
    easyloggingpp
    GIT_REPOSITORY https://github.com/muflihun/easyloggingpp.git
    GIT_TAG        v9.97.0
)
FetchContent_MakeAvailable(easyloggingpp)

# Include easylogging++ headers
include_directories(${easyloggingpp_SOURCE_DIR}/src)

# Compile `easylogging++.cc` as a part of the build
add_library(easyloggingpp STATIC ${easyloggingpp_SOURCE_DIR}/src/easylogging++.cc)

# Add executables from all .cxx files in scripts
file(GLOB scripts ${CMAKE_CURRENT_SOURCE_DIR}/scripts/*.cxx)
foreach(script ${scripts})
    get_filename_component(script_name ${script} NAME_WE)
    add_executable(${script_name} ${script})

    # Ensure easylogging++ headers are available
    target_include_directories(${script_name} PUBLIC ${easyloggingpp_SOURCE_DIR}/src)

    # Link easylogging++ (compiled) and argparse
    target_link_libraries(${script_name} easyloggingpp argparse ${ROOT_LIBRARIES})
endforeach(script)
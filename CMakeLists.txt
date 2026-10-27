cmake_minimum_required(VERSION 3.16)
project(kawasaki_model LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(kawasaki_core INTERFACE)
target_include_directories(kawasaki_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(kawasaki src/kawasaki.cpp)
target_link_libraries(kawasaki PRIVATE kawasaki_core)

if(MSVC)
    target_compile_options(kawasaki PRIVATE /W4)
else()
    target_compile_options(kawasaki PRIVATE -Wall -Wextra -Wpedantic -O2)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(test_lattice tests/test_lattice.cpp)
    target_link_libraries(test_lattice PRIVATE kawasaki_core)
    add_test(NAME lattice_energy_and_magnetization COMMAND test_lattice)
    add_test(NAME simulator_help COMMAND kawasaki --help)
endif()

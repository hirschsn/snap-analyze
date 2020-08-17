cmake_minimum_required(VERSION 3.10)

project(snap-analyze)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(DEFAULT_COMPILE_OPTIONS "-Wall" "-Wextra" "-pedantic")

add_executable(sana-bond sana-bond.cpp)
target_compile_options(sana-bond PRIVATE ${DEFAULT_COMPILE_OPTIONS})

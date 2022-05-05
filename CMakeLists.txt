cmake_minimum_required(VERSION 3.21)
project(ViralSpread)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ViralSpread main.cpp)
target_compile_features(ViralSpread PRIVATE cxx_std_17)
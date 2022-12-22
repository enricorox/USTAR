cmake_minimum_required(VERSION 3.24)
project(USTAR)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-O3 -W)

add_executable(ustar src/main.cpp src/DBG.cpp src/SPSS.cpp)

add_executable(ustar-test src/ustar-test.cpp src/DBG.cpp src/SPSS.cpp)

cmake_minimum_required(VERSION 3.16)
project(USTAR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)

#add_compile_options(-O0 -W)
add_compile_options(-O3 -W)

add_executable(ustar
        src/ustar.cpp
        src/DBG.cpp
        src/SPSS.cpp
        src/Sorter.cpp
        src/Encoder.cpp
        src/commons.cpp
        src/bwt.hpp
        )

add_executable(ustar-test
        src/ustar-test.cpp
        src/DBG.cpp
        src/SPSS.cpp
        src/Sorter.cpp
        src/Encoder.cpp
        src/commons.cpp
        src/bwt.hpp
        )

add_executable(ustarx
        src/ustarx.cpp
        src/Decoder.cpp
        src/DBG.cpp
        src/commons.cpp
        src/bwt.hpp
        )

add_executable(ustars
        src/ustars.cpp
        src/Analyzer.cpp
        )
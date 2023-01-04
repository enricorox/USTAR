cmake_minimum_required(VERSION 3.16)
project(USTAR)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-O3 -W)

add_executable(ustar
        src/ustar.cpp
        src/DBG.cpp src/DBG.h
        src/SPSS.cpp src/SPSS.h
        src/Sorter.cpp src/Sorter.h
        src/Encoder.cpp src/Encoder.h
        src/commons.cpp src/commons.h
        src/bwt.hpp
        src/Postprocessor.cpp src/Postprocessor.h
        )

add_executable(ustar-test src/ustar-test.cpp src/DBG.cpp src/SPSS.cpp src/Sorter.cpp src/Encoder.cpp src/commons.cpp src/bwt.hpp src/Postprocessor.cpp)

add_executable(ustarx src/ustarx.cpp src/Decoder.cpp src/DBG.cpp src/commons.cpp src/bwt.hpp)
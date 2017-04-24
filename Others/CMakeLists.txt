cmake_minimum_required(VERSION 3.7)
project(Others)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        computeMOWSEMatrix.cpp
        msmsDBSearch.cpp
        pmfDBSearch.cpp)

add_executable(Others ${SOURCE_FILES})
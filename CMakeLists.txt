cmake_minimum_required(VERSION 2.8.9)

project(test)
set(SOURCES test.cpp IMParton.cpp)

add_library(IMParton SHARED IMParton.cpp)
add_executable(testIMParton ${SOURCES})

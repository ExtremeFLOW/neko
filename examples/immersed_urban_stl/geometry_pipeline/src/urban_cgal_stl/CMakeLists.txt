cmake_minimum_required(VERSION 3.16)
project(urban_cgal_stl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)

add_executable(urban_cgal_stl urban_cgal_stl.cpp)
target_link_libraries(urban_cgal_stl PRIVATE CGAL::CGAL)

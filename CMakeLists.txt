cmake_minimum_required(VERSION 3.10)
project(Solving_Systems_Of_Linear_Equations)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(Solving_Systems_Of_Linear_Equations
    main.cpp
    src/matrix.cpp
    src/linear_system.cpp
    src/utils.cpp
)
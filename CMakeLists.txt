cmake_minimum_required(VERSION 3.12)
project(CMPT431)

set(CMAKE_CXX_STANDARD 11)



add_executable(ass0
        assignment0/driver.cpp
#        assignment0/solution.cpp
        )


target_link_libraries(ass0 pthread)

add_executable(ass1
        assignment1/pi_calculation.cpp
        #        assignment0/solution.cpp
        )


target_link_libraries(ass1 pthread)


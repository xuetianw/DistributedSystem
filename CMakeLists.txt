cmake_minimum_required(VERSION 3.12)
project(CMPT431)

set(CMAKE_CXX_STANDARD 11)



add_executable(ass0
        assignment0/driver.cpp
#        assignment0/solution.cpp
        )


target_link_libraries(ass0 pthread)

add_executable(ass1_part1
        assignment1/pi_calculation.cpp
        )

target_link_libraries(ass1_part1 pthread)

add_executable(pi_calculation_parallel
        assignment1/pi_calculation_parallel.cpp
        )

target_link_libraries(pi_calculation_parallel pthread)

add_executable(triangle_counting_parallel
        assignment1/triangle_counting_parallel.cpp
        )

target_link_libraries(triangle_counting_parallel pthread)

add_executable(page_rank_parallel
        assignment2/page_rank_parallel.cpp
        #        assignment0/solution.cpp
        )


target_link_libraries(page_rank_parallel pthread)


add_executable(test
        assignment1/test.cpp
        #        assignment0/solution.cpp
        )


target_link_libraries(test pthread)


add_executable(test2
        assignment2/test.cpp
        #        assignment0/solution.cpp
        )


target_link_libraries(test2 pthread)

add_executable(AtomicTest
        assignment2/AtomicTest.cpp
        #        assignment0/solution.cpp
        )


target_link_libraries(AtomicTest pthread)



add_executable(page_rank_parallel_atomic
        assignment2/page_rank_parallel_atomic.cpp
        #        assignment0/solution.cpp
        )

target_link_libraries(page_rank_parallel_atomic pthread)



add_executable(triangle_counting_parallel_ass3
        assignment3/triangle_counting_parallel.cpp
        #        assignment0/solution.cpp
        )

target_link_libraries(triangle_counting_parallel_ass3 pthread)


add_executable(page_rank_parallel_atomic_ass3
        assignment3/page_rank.cpp
        #        assignment0/solution.cpp
        )

target_link_libraries(page_rank_parallel_atomic_ass3 pthread)


add_executable(triangle_counting_parallel_ass4
        assignment4/triangle_counting_parallel_ass4.cpp
        #        assignment0/solution.cpp
        )

target_link_libraries(triangle_counting_parallel_ass4 pthread)

add_executable(page_rank_parallel_ass4
        assignment4/page_rank_parallel_ass4.cpp
        #        assignment0/solution.cpp
        )

target_link_libraries(page_rank_parallel_ass4 pthread)
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "12345678"

uint c_const = (uint) RAND_MAX + (uint) 1;

std::mutex mtx;           // mutex for critical section

inline double get_random_coordinate(uint* random_seed) {
    return ((double) rand_r(random_seed)) / c_const;
}

struct res_data {
    uint* points_to_generate_arr;
    uint* circle_points;
    double* time_taken_s;
    uint total_points_in_circle;
};

struct arg_struct {
    arg_struct(res_data& resData) : resData(resData) {}

    res_data& resData;
};


void* calculatingFunction(void* data, int thread_id) {
    timer serial_timer;
    serial_timer.start();

    arg_struct& args = *(arg_struct*) data;;
    res_data& resData = args.resData;
    uint points_to_generate = resData.points_to_generate_arr[thread_id];
    uint random_seed = thread_id;

    uint circle_count = 0;
    double x_coord, y_coord;
    for (uint i = 0; i < points_to_generate; i++) {
        x_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
        y_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
        if ((sqr(x_coord) + sqr(y_coord)) <= 1.0) {
            circle_count++;
        }
    }


    mtx.lock();
    resData.circle_points[thread_id] = circle_count;
    resData.total_points_in_circle += circle_count;
    double time_taken = serial_timer.stop();
    resData.time_taken_s[thread_id] = time_taken;

    mtx.unlock();
}


void piCalculation(uint n, uint n_workers) {
    timer serial_timer;
    double time_taken = 0.0;

    serial_timer.start();
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------

//    std::thread threads[n_workers];

    res_data resData{};
    resData.points_to_generate_arr = new uint[n_workers];

    for (int i = 0; i < n_workers - 1; i++) {
        resData.points_to_generate_arr[i] = n / n_workers;
    }
    resData.points_to_generate_arr[n_workers - 1] =
            n % n_workers == 0 ? n / n_workers : n - (n / n_workers) * (n_workers - 1);

    resData.circle_points = new uint[n_workers];
    resData.time_taken_s = new double[n_workers];

    arg_struct argS{resData};
    std::vector<std::thread> threads;
    threads.reserve(n_workers);
    for (uint i = 0; i < n_workers; i++)
        threads.emplace_back(calculatingFunction, &argS, i);

    // -------------------------------------------------------------------

    std::cout << "thread_id, points_generated, circle_points, time_taken\n";
    // Print the above statistics for each thread
    // Example output for 2 threads:
    // thread_id, points_generated, circle_points, time_taken
    // 1, 100, 90, 0.12
    // 0, 100, 89, 0.12

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
//        pthread_join(pthread_ts[i], nullptr);
        printf("%d,\t %d,\t %d,\t %f\n", i, resData.points_to_generate_arr[i], resData.circle_points[i],
               resData.time_taken_s[i]);
    }

    double pi_value = 4.0 * (double) resData.total_points_in_circle / (double) n;

    time_taken = serial_timer.stop();
    // Print the overall statistics
    std::cout << "Total points generated : " << n << "\n";
    std::cout << "Total points in circle : " << resData.total_points_in_circle << "\n";
    std::cout << "Result : " << std::setprecision(VAL_PRECISION) << pi_value
              << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";

    delete[] resData.circle_points;
    delete[] resData.time_taken_s;
    delete[] resData.points_to_generate_arr;
}

int main(int argc, char* argv[]) {


    // Initialize command line arguments
    cxxopts::Options options("pi_calculation",
                             "Calculate pi using serial and parallel execution");
    options.add_options(
            "custom",
            {
                    {"nPoints",  "Number of points",
                            cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
                    {"nWorkers", "Number of workers",
                            cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_points = cl_options["nPoints"].as<uint>();
    uint n_workers = cl_options["nWorkers"].as<uint>();
//    uint n_points = 123456001;
//    uint n_workers = 100;
    std::cout << std::fixed;
    std::cout << "Number of points : " << n_points << "\n";
    std::cout << "Number of workers : " << n_workers << "\n";

    piCalculation(n_points, n_workers);

    return 0;
}

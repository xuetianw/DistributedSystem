#include "core/graph.h"
#include "core/utils.h"
#include <future>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

std::mutex mtx;           // mutex for critical section

uintV countTriangles(uintV* array1, uintE len1, uintV* array2, uintE len2,
                     uintV u, uintV v) {

    uintE i = 0, j = 0; // indexes for array1 and array2
    uintV count = 0;

    if (u == v)
        return count;

    while ((i < len1) && (j < len2)) {
        if (array1[i] == array2[j]) {
            if ((array1[i] != u) && (array1[i] != v)) {
                count++;
            }
            i++;
            j++;
        } else if (array1[i] < array2[j]) {
            i++;
        } else {
            j++;
        }
    }
    return count;
}


struct res_data {
    long* triangle_count_arr;
    long* vertices_to_calculate_arr;
    double* time_taken_s_arr;
    long triangles_total;
};

struct arg_struct {
    arg_struct(res_data& resData, Graph& graph) : graph(graph), resData(resData) {}

    Graph& graph;
    res_data& resData;
};


void* calculatingFunction(void* data, int thread_id) {
    timer serial_timer;
    serial_timer.start();

    arg_struct& args = *(arg_struct*) data;;
    Graph& g = args.graph;
    res_data& resData = args.resData;

    uint num_of_vertices = resData.vertices_to_calculate_arr[thread_id];

    uintV starting_vertex = thread_id == 0 ? 0 :
            thread_id * resData.vertices_to_calculate_arr[thread_id - 1];

    uintV finishing_vex = starting_vertex + num_of_vertices - 1;
    uint triangle_count = 0;

//
    for (uintV u = starting_vertex; u <= finishing_vex; u++) {
        // For each outNeighbor v, find the intersection of inNeighbor(u) and
        // outNeighbor(v)
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(), u, v);
        }
    }


    mtx.lock();

    resData.triangle_count_arr[thread_id] = triangle_count;
    resData.triangles_total += triangle_count;

    double time_taken = serial_timer.stop();
    resData.time_taken_s_arr[thread_id] = time_taken;

    mtx.unlock();

    delete (arg_struct*) data;
}


void triangleCountSerial(Graph& g, uint n_workers) {
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;
    t1.start();

    // The outNghs and inNghs for a given vertex are already sorted

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------

    res_data resData{};
    resData.time_taken_s_arr = new double[n_workers];
    resData.vertices_to_calculate_arr = new long[n_workers];
    resData.triangle_count_arr = new long[n_workers];

    for (int i = 0; i < n_workers - 1; i++) {
        resData.vertices_to_calculate_arr[i] = n / n_workers;
    }
    resData.vertices_to_calculate_arr[n_workers - 1] =
            n % n_workers == 0 ? n / n_workers : n - (n / n_workers) * (n_workers - 1);


    std::thread threads[n_workers];

    arg_struct argS{resData, g};
    for (int i = 0; i < n_workers; i++) {
        threads[i] = std::thread{calculatingFunction, &argS, i};
    }


    std::cout << "thread_id, triangle_count, time_taken\n";

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
//        pthread_join(pthread_ts[i], nullptr);
        printf("%d,\t %ld,\t %f\n", i, resData.triangle_count_arr[i], resData.time_taken_s_arr[i]);
    }
    // -------------------------------------------------------------------
    // Here, you can just print the number of non-unique triangles counted by each
    // thread std::cout << "thread_id, triangle_count, time_taken\n"; Print the
    // above statistics for each thread Example output for 2 threads: thread_id,
    // triangle_count, time_taken 1, 102, 0.12 0, 100, 0.12

    time_taken = t1.stop();
    // Print the overall statistics
    std::cout << "Number of triangles : " << resData.triangles_total << "\n";
    std::cout << "Number of unique triangles : " << resData.triangles_total / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";

    delete[] resData.triangle_count_arr;
    delete[] resData.time_taken_s_arr;
    delete[] resData.vertices_to_calculate_arr;
}

int main(int argc, char* argv[]) {
    cxxopts::Options options(
            "triangle_counting_serial",
            "Count the number of triangles using serial and parallel execution");
    options.add_options(
            "custom",
            {
                    {"nWorkers",  "Number of workers",
                            cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                    {"inputFile", "Input graph file path",
                            cxxopts::value<std::string>()->default_value(
                                    "/scratch/assignment1/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
//    std::string input_file_path = cl_options["inputFile"].as<std::string>();
//    uint n_workers = 10;
    std::string input_file_path = "/home/fred/DistributedSystemProject/assignment1/inputfiles/lj";
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    triangleCountSerial(g, n_workers);

    return 0;
}

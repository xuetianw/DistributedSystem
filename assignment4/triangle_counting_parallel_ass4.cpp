#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>


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
            } else {
                // triangle with self-referential edge -> ignore
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


std::mutex mtx;           // mutex for critical section
std::mutex buffer_mtx;

struct res_data {
    long* triangle_count_arr;
    long* edge_arr;
    long* num_vertices_arr;
    double* time_taken_s_arr;
    long triangles_total;
    std::atomic<int32_t> starting_index;

    int32_t getNextVertexToBeProcessed() {
        buffer_mtx.lock();
        int32_t copy = --starting_index;
        buffer_mtx.unlock();
        return copy;
    }
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

    uint triangle_count = 0;


    while (true) {
        uintE u = resData.getNextVertexToBeProcessed();
        if (u <= -1) break;
        resData.num_vertices_arr[thread_id]++;
        uintE out_degree = g.vertices_[u].getOutDegree();

        resData.edge_arr[thread_id] += out_degree;

        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(), u, v);
        }
    }

    resData.triangle_count_arr[thread_id] = triangle_count;

    mtx.lock();

    resData.triangles_total += triangle_count;

    mtx.unlock();

    double time_taken = serial_timer.stop();
    resData.time_taken_s_arr[thread_id] = time_taken;
}

void triangleVertexBasedCountParallel(Graph& g, uint n_workers) {
    timer t1;
    t1.start();

    res_data resData{};
    resData.time_taken_s_arr = new double[n_workers];
    resData.num_vertices_arr = new long[n_workers];
    resData.triangle_count_arr = new long[n_workers];
    resData.edge_arr = new long[n_workers];

    for (int i = 0; i < n_workers; i++) {
        resData.time_taken_s_arr[i] = 0;
        resData.num_vertices_arr[i] = 0;
        resData.triangle_count_arr[i] = 0;
        resData.edge_arr[i] = 0;
    }

    arg_struct argS{resData, g};

    std::vector<std::thread> threads;
    threads.reserve(n_workers);

    resData.starting_index = g.n_;

    for (int i = 0; i < n_workers; i++)
        threads.emplace_back(calculatingFunction, &argS, i);

    std::cout << "thread_id, num_vertices, num_edges, triangle_count, time_taken\n";

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
        printf("%d,\t %ld,\t %ld,\t %ld,\t %f\n", i, resData.num_vertices_arr[i],
               resData.edge_arr[i], resData.triangle_count_arr[i],
               resData.time_taken_s_arr[i]);
    }

    double time_taken = t1.stop();

    std::cout << "Number of triangles : " << resData.triangles_total << "\n";
    std::cout << "Number of unique triangles : " << resData.triangles_total / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
    delete[] resData.triangle_count_arr;
    delete[] resData.time_taken_s_arr;
    delete[] resData.edge_arr;
    delete[] resData.num_vertices_arr;
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
                                    "/scratch/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    triangleVertexBasedCountParallel(g, n_workers);
    return 0;
}

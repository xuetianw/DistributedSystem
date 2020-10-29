#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <set>

std::mutex mtx;           // mutex for critical section
std::mutex mtx2;


uintV countTriangles(uintV* array1, uintE len1, uintV* array2, uintE len2,
                     uintV u, uintV v);

// strategy 2 edge-based

struct res_data {
    long* triangle_count_arr;
    long* edge_arr;
    long* num_vertices_arr;
    double* time_taken_s_arr;
    long triangles_total;
};


struct edge_arg_struct {
    edge_arg_struct(res_data& resData, Graph& graph) : graph(graph), resData(resData) {}

    Graph& graph;
    res_data& resData;
};


void* edgeCalculatingFunction(void* data, int thread_id, std::vector<std::pair<uintV, uintV>>* edges) {
    timer serial_timer;
    serial_timer.start();


    edge_arg_struct& args = *(edge_arg_struct*) data;;
    Graph& g = args.graph;
    res_data& resData = args.resData;

    uint num_of_vertices = resData.edge_arr[thread_id];

    uintV starting_edge = thread_id == 0 ? 0 :
                          thread_id * resData.edge_arr[thread_id - 1];

    uintV finishing_edge = starting_edge + num_of_vertices - 1;
    uint triangle_count = 0;

    for (uintV u = starting_edge; u <= finishing_edge; u++) {
        uintV vertex1 = edges->at(u).first;
        uintV vertex2 = edges->at(u).second;

        triangle_count += countTriangles(g.vertices_[vertex1].getInNeighbors(),
                                         g.vertices_[vertex1].getInDegree(),
                                         g.vertices_[vertex2].getOutNeighbors(),
                                         g.vertices_[vertex2].getOutDegree(), vertex1, vertex2);
    }

    mtx.lock();

    resData.triangle_count_arr[thread_id] = triangle_count;
    resData.triangles_total += triangle_count;

    double time_taken = serial_timer.stop();
    resData.time_taken_s_arr[thread_id] = time_taken;
    resData.num_vertices_arr[thread_id] = 0;

    mtx.unlock();
}

void triangleCountEdgeBasedParallel(Graph& g, uint n_workers) {
    std::vector<std::pair<uintV, uintV>> edges;

    uintV n = g.n_;
    for (uintV u = 0; u < n; u++) {
        Vertex& vertex = g.vertices_[u];
        uintE out_degree = vertex.getOutDegree();
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = vertex.getOutNeighbor(i);
            edges.emplace_back(u, v);
        }
    }

    timer t1;
    t1.start();

    uintV edge_count = edges.size();

    res_data edgeResData{};
    edgeResData.time_taken_s_arr = new double[n_workers];
    edgeResData.edge_arr = new long[n_workers];
    edgeResData.num_vertices_arr = new long[n_workers];
    edgeResData.triangle_count_arr = new long[n_workers];

    timer t2;
    double partitionTime = 0.0;

    for (int i = 0; i < n_workers - 1; i++) {
        t2.start();
        edgeResData.edge_arr[i] = edge_count / n_workers;
        partitionTime += t2.stop();
    }

    t2.start();

    edgeResData.edge_arr[n_workers - 1] =
            edge_count % n_workers == 0 ? edge_count / n_workers : edge_count -
                                                                   (edge_count / n_workers) * (n_workers - 1);

    partitionTime += t2.stop();

    edge_arg_struct argS{edgeResData, g};
    std::thread threads[n_workers];
    for (int i = 0; i < n_workers; i++)
        threads[i] = std::thread{edgeCalculatingFunction, &argS, i, &edges};

    std::cout << "thread_id, num_vertices, num_edges, triangle_count, time_taken\n";

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
        printf("%d,\t %ld,\t %ld,\t %ld,\t %f\n", i, edgeResData.num_vertices_arr[i],
               edgeResData.edge_arr[i], edgeResData.triangle_count_arr[i],
               edgeResData.time_taken_s_arr[i]);
    }

    double time_taken = t1.stop();
    // Print the overall statistics
    std::cout << "Number of triangles : " << edgeResData.triangles_total << "\n";
    std::cout << "Number of unique triangles : " << edgeResData.triangles_total / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
    std::cout << "Partitioning time (in seconds) : " << partitionTime << "\n";
    delete[] edgeResData.triangle_count_arr;
    delete[] edgeResData.time_taken_s_arr;
    delete[] edgeResData.edge_arr;
    delete[] edgeResData.num_vertices_arr;
}


// strategy 1


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

    uint num_of_vertices = resData.num_vertices_arr[thread_id];

    uintV starting_vertex = thread_id == 0 ? 0 :
                            thread_id * resData.num_vertices_arr[thread_id - 1];

    uintV finishing_vex = starting_vertex + num_of_vertices - 1;
    uint triangle_count = 0;

    for (uintV u = starting_vertex; u <= finishing_vex; u++) {
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

    resData.triangles_total = resData.triangles_total + 1;
    mtx.unlock();

    double time_taken = serial_timer.stop();
    resData.time_taken_s_arr[thread_id] = time_taken;

}

// strategy 1 vertices-based

void triangleVertexBasedCountParallel(Graph& g, uint n_workers) {
    uintV n = g.n_;
    timer t1;
    t1.start();

    res_data resData{};
    resData.time_taken_s_arr = new double[n_workers];
    resData.num_vertices_arr = new long[n_workers];
    resData.triangle_count_arr = new long[n_workers];
    resData.edge_arr = new long[n_workers];

    timer t2;
    double partitionTime = 0.0;
    for (int i = 0; i < n_workers - 1; i++) {
        t2.start();
        resData.num_vertices_arr[i] = n / n_workers;
        partitionTime += t2.stop();
    }
    t2.start();
    resData.num_vertices_arr[n_workers - 1] =
            n % n_workers == 0 ? n / n_workers : n - (n / n_workers) * (n_workers - 1);
    partitionTime += t2.stop();

    arg_struct argS{resData, g};

    std::vector<std::thread> threads;
    threads.reserve(n_workers);
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
    // Print the overall statistics
    std::cout << "Number of triangles : " << resData.triangles_total << "\n";
    std::cout << "Number of unique triangles : " << resData.triangles_total / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
    std::cout << "Partitioning time (in seconds) : " << partitionTime << "\n";
    delete[] resData.triangle_count_arr;
    delete[] resData.time_taken_s_arr;
    delete[] resData.edge_arr;
    delete[] resData.num_vertices_arr;
}


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

void triangleCountSerial(Graph& g) {
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken = 0.0;
    timer t1;
    t1.start();
    for (uintV u = 0; u < n; u++) {
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(), u, v);
        }
    }
    time_taken = t1.stop();
    std::cout << "Number of triangles : " << triangle_count << "\n";
    std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
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
                    {"strategy",  "Strategy to be used",
                            cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                    {"inputFile", "Input graph file path",
                            cxxopts::value<std::string>()->default_value(
                                    "/scratch/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    switch (strategy) {
        case 0:
//            std::cout << "\nSerial\n";
            triangleCountSerial(g);
            break;
        case 1:
//            std::cout << "\nVertex-based work partitioning\n";
            triangleVertexBasedCountParallel(g, n_workers);
            break;
        case 2:
//            std::cout << "\nEdge-based work partitioning\n";
            triangleCountEdgeBasedParallel(g, n_workers);
            break;
        default:
            break;
    }

    return 0;
}


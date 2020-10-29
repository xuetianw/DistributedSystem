#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <set>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

// strategy 1 Parallel vertices-based

struct res_data {
    uintV* num_vertices_arr;
    uintV* edge_arr;
    uintV* edges_processed;
    uintV* barrier1_time_arr;
    uintV* barrier2_time_arr;
    double* time_taken_s;
    uintV* vertices_processed;
    uintV max_iterations;
};

struct arg_struct {
    arg_struct(res_data& resData, CustomBarrier& b, Graph& graph) : resData(resData), graph(graph), b(b) {}

    Graph& graph;
    CustomBarrier& b;
    res_data& resData;
};

std::mutex mutex;

void*
calculatingFunction(std::atomic<PageRankType>* pr_curr, std::atomic<PageRankType>* pr_next, void* data, int thread_id) {

    timer serial_timer;
    serial_timer.start();

    arg_struct& args = *(arg_struct*) data;;
    Graph& g = args.graph;
    res_data& resData = args.resData;
    uintV max_iterations = resData.max_iterations;

    uint num_of_vertices = resData.num_vertices_arr[thread_id];

    uintV starting_vertex = thread_id == 0 ? 0 :
                            thread_id * resData.num_vertices_arr[thread_id - 1];
    uintV ending_vertex = starting_vertex + num_of_vertices - 1;

    CustomBarrier& b = args.b;

    timer t;
    t.start();
    double barrier1_time = 0;
    double barrier2_time = 0;

    resData.edge_arr[thread_id] = 0;
    resData.vertices_processed[thread_id] = 0; // initialized to 0, c++ does not have default 0 value

    for (int i = 0; i < max_iterations; i++) {

        resData.vertices_processed[thread_id] += resData.num_vertices_arr[thread_id];

        for (uintV u = starting_vertex; u <= ending_vertex; u++) {
            uintE out_degree = g.vertices_[u].getOutDegree();

            resData.edge_arr[thread_id] += out_degree;

            for (uintE j = 0; j < out_degree; j++) {
                uintV v = g.vertices_[u].getOutNeighbor(j);

                float temp = pr_next[v];
                while (!pr_next[v].compare_exchange_weak(temp, pr_next[v] + pr_curr[u] / out_degree));
            }
        }
        b.wait();
        barrier1_time += t.stop();
        t.start();

        for (uintV v = starting_vertex; v <= ending_vertex; v++) {

//            resData.vertices_processed[thread_id]++;

            pr_next[v] = PAGE_RANK(pr_next[v]);
            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v].load();
            pr_next[v] = 0.0;
        }

        b.wait();
        barrier2_time += t.stop();
    }

    resData.barrier2_time_arr[thread_id] = barrier2_time;
    resData.barrier1_time_arr[thread_id] = barrier1_time;

    double time_taken = serial_timer.stop();
    resData.time_taken_s[thread_id] = time_taken;
}


void pageRankVertexBasedParallel(Graph& g, int max_iterations, int n_workers) {
    uintV n = g.n_;
    double time_taken = 0.0;
    timer t1;
    t1.start();

    std::atomic<PageRankType>* pr_curr = new std::atomic<PageRankType>[n];
    std::atomic<PageRankType>* pr_next = new std::atomic<PageRankType>[n];

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    std::thread threads[n_workers];
    CustomBarrier b{n_workers};

    res_data resData{};
    resData.max_iterations = max_iterations;
    resData.num_vertices_arr = new uintV[n_workers];
    resData.vertices_processed = new uintV[n_workers];
    resData.time_taken_s = new double[n_workers];
    resData.edge_arr = new int32_t[n_workers];
    resData.barrier1_time_arr = new int32_t[n_workers];
    resData.barrier2_time_arr = new int32_t[n_workers];

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

    arg_struct aStruct{resData, b, g};

    for (int i = 0; i < n_workers; i++)
        threads[i] = std::thread{calculatingFunction, pr_curr, pr_next, &aStruct, i};

    std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, total_time\n";

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
        printf("%d,\t %d,\t %d,\t %d,\t %d,\t %f\n", i, resData.vertices_processed[i],
               resData.edge_arr[i], resData.barrier1_time_arr[i], resData.barrier2_time_arr[i],
               resData.time_taken_s[i]);
    }


    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
//        std::cerr << sum_of_page_ranks << std::endl;
    }
    time_taken = t1.stop();
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Partitioning time (in seconds) : " << partitionTime << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    delete[] pr_curr;
    delete[] pr_next;

    delete[] resData.vertices_processed;
    delete[] resData.time_taken_s;
    delete[] resData.num_vertices_arr;
    delete[] resData.edge_arr;
    delete[] resData.barrier1_time_arr;
    delete[] resData.barrier2_time_arr;
}



// strategy 2 Edge-Based -based

void*
edgeCalculatingFunction(std::atomic<PageRankType>* pr_curr, std::atomic<PageRankType>* pr_next, void* data,
                        int thread_id, std::vector<std::pair<uintV, uintV>>* edges) {

    timer serial_timer;
    serial_timer.start();

    arg_struct& args = *(arg_struct*) data;;
    Graph& g = args.graph;
    res_data& resData = args.resData;
    uintV max_iterations = resData.max_iterations;

    uint num_of_edges = resData.edge_arr[thread_id];
    uintV starting_edge = thread_id == 0 ? 0 :
                          thread_id * resData.edge_arr[thread_id - 1];
    uintV finishing_edge = starting_edge + num_of_edges - 1;

    CustomBarrier& b = args.b;

    timer t;
    t.start();
    double barrier1_time = 0;
    double barrier2_time = 0;

    resData.num_vertices_arr[thread_id] = 0;
    resData.edges_processed[thread_id] = 0;

    for (int i = 0; i < max_iterations; i++) {

        resData.edges_processed[thread_id] += resData.edge_arr[thread_id];

        for (uintV u = starting_edge; u <= finishing_edge; u++) {

            uintV vertex1 = edges->at(u).first;
            uintV vertex2 = edges->at(u).second;

            float temp = pr_next[vertex2];
            uintE out_degree = g.vertices_[vertex1].getOutDegree();
            while (!pr_next[vertex2].compare_exchange_weak(temp, pr_next[vertex2] + pr_curr[vertex1] / out_degree));

        }

        b.wait();

        barrier1_time += t.stop();
        t.start();


        if (thread_id == 0) {
            for (uintV v = 0; v < g.n_; v++) {

                resData.num_vertices_arr[0]++;

                pr_next[v] = PAGE_RANK(pr_next[v]);
                // reset pr_curr for the next iteration
                pr_curr[v] = pr_next[v].load();
                pr_next[v] = 0.0;
            }
        }

        b.wait();
        barrier2_time += t.stop();
    }


    resData.barrier2_time_arr[thread_id] = barrier2_time;
    resData.barrier1_time_arr[thread_id] = barrier1_time;

    double time_taken = serial_timer.stop();
    resData.time_taken_s[thread_id] = time_taken;

}


void pageRankEdgeBasedParallel(Graph& g, int max_iterations, int n_workers) {
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

    std::atomic<PageRankType>* pr_curr = new std::atomic<PageRankType>[n];
    std::atomic<PageRankType>* pr_next = new std::atomic<PageRankType>[n];


    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }
    std::thread threads[n_workers];
    CustomBarrier b{n_workers};

    res_data resData{};
    resData.max_iterations = max_iterations;
    resData.num_vertices_arr = new uintV[n_workers];
    resData.time_taken_s = new double[n_workers];
    resData.edge_arr = new int32_t[n_workers];
    resData.barrier1_time_arr = new int32_t[n_workers];
    resData.barrier2_time_arr = new int32_t[n_workers];
    resData.edges_processed = new int32_t[n_workers];

    timer t2;
    double partitionTime = 0.0;
    for (int i = 0; i < n_workers - 1; i++) {
        t2.start();
        resData.edge_arr[i] = edge_count / n_workers;
        partitionTime += t2.stop();
    }
    t2.start();

    resData.edge_arr[n_workers - 1] =
            edge_count % n_workers == 0 ? edge_count / n_workers : edge_count -
                                                                   (edge_count / n_workers) * (n_workers - 1);
    partitionTime += t2.stop();

    arg_struct aStruct{resData, b, g};

    for (int i = 0; i < n_workers; i++)
        threads[i] = std::thread{edgeCalculatingFunction, pr_curr, pr_next, &aStruct, i, &edges};

    std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, total_time\n";

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
        printf("%d,\t %d,\t %d,\t %d,\t %d,\t %f\n", i, resData.num_vertices_arr[i],
               resData.edges_processed[i], resData.barrier1_time_arr[i], resData.barrier2_time_arr[i],
               resData.time_taken_s[i]);
    }


    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
//        std::cerr << sum_of_page_ranks << std::endl;
    }
    double time_taken = t1.stop();
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Partitioning time (in seconds) : " << partitionTime << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    delete[] pr_curr;
    delete[] pr_next;

    delete[] resData.time_taken_s;
    delete[] resData.num_vertices_arr;
    delete[] resData.edge_arr;
    delete[] resData.barrier1_time_arr;
    delete[] resData.barrier2_time_arr;
    delete[] resData.edges_processed;
}



// strategy 0 Serial

void pageRankSerial(Graph& g, int max_iters) {
    uintV n = g.n_;

    PageRankType* pr_curr = new PageRankType[n];
    PageRankType* pr_next = new PageRankType[n];

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    timer t1;
    double time_taken = 0.0;
    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int iter = 0; iter < max_iters; iter++) {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = 0; u < n; u++) {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }
        for (uintV v = 0; v < n; v++) {
            pr_next[v] = PAGE_RANK(pr_next[v]);
            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";
    delete[] pr_curr;
    delete[] pr_next;
}

int main(int argc, char* argv[]) {
    cxxopts::Options options(
            "page_rank_push",
            "Calculate page_rank using serial and parallel execution");
    options.add_options(
            "",
            {
                    {"nWorkers",    "Number of workers",
                            cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
                    {"nIterations", "Maximum number of iterations",
                            cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                    {"strategy",    "Strategy to be used",
                            cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                    {"inputFile",   "Input graph file path",
                            cxxopts::value<std::string>()->default_value(
                                    "/scratch/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint strategy = cl_options["strategy"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Task decomposition strategy : " << strategy << "\n";
    std::cout << "Iterations : " << max_iterations << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";
    switch (strategy) {
        case 0:
//            std::cout << "\nSerial\n";
            pageRankSerial(g, max_iterations);
            break;
        case 1:
//            std::cout << "\nVertex-based work partitioning\n";
            pageRankVertexBasedParallel(g, max_iterations, n_workers);
            break;
        case 2:
//            std::cout << "\nEdge-based work partitioning\n";
            pageRankEdgeBasedParallel(g, max_iterations, n_workers);
            break;
        default:
            break;
    }

    return 0;
}

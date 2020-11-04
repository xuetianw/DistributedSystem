#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

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

std::mutex buffer_mtx;

struct res_data {
    uintV* edge_arr;
    uintV* edges_processed;
    uintV* barrier1_time_arr;
    uintV* barrier2_time_arr;
    double* time_taken_s;
    uintV* vertices_processed;
    uintV max_iterations;

    int32_t starting_index;

    int32_t getNextVertexToBeProcessed() {
        buffer_mtx.lock();
        int32_t copy = --starting_index;
        buffer_mtx.unlock();
        return copy;
    }
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

    CustomBarrier& b = args.b;

    timer t;
    t.start();
    double barrier1_time = 0;
    double barrier2_time = 0;

    resData.edge_arr[thread_id] = 0;
    resData.vertices_processed[thread_id] = 0; // initialized to 0, c++ does not have default 0 value

    for (int i = 0; i < max_iterations; i++) {

        while (true) {
            uintE u = resData.getNextVertexToBeProcessed();
            if (u <= -1) break;

            resData.vertices_processed[thread_id]++;
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
        if (thread_id == 0)
            resData.starting_index = g.n_;

        b.wait();
        t.start();

        while (true) {
            uintE v = resData.getNextVertexToBeProcessed();
            if (v <= -1) break;

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


void pageRankParallel(Graph& g, uint max_iterations, int n_workers) {
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
    resData.vertices_processed = new uintV[n_workers];
    resData.time_taken_s = new double[n_workers];
    resData.edge_arr = new int32_t[n_workers];
    resData.barrier1_time_arr = new int32_t[n_workers];
    resData.barrier2_time_arr = new int32_t[n_workers];

    resData.starting_index = g.n_;

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
    for (uintV u = 0; u < n; u++)
            sum_of_page_ranks += pr_curr[u];
    time_taken = t1.stop();
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    delete[] pr_curr;
    delete[] pr_next;

    delete[] resData.vertices_processed;
    delete[] resData.time_taken_s;
    delete[] resData.edge_arr;
    delete[] resData.barrier1_time_arr;
    delete[] resData.barrier2_time_arr;

}


//  Granularity

struct res_data_g {
    uintV* edge_arr;
    uintV* edges_processed;
    uintV* barrier1_time_arr;
    uintV* barrier2_time_arr;
    double* time_taken_s;
    uintV* vertices_processed;
    uintV max_iterations;

    uintV granularity;

    int32_t starting_index;

    int32_t getNextVertexToBeProcessed() {
        buffer_mtx.lock();
        int copy = starting_index;
        starting_index -= granularity;
        buffer_mtx.unlock();
        return copy;
    }
};


struct arg_struct_g {
    arg_struct_g(res_data_g& resData, CustomBarrier& b, Graph& graph) : graph(graph), b(b), resData(resData) {}

    Graph& graph;
    CustomBarrier& b;
    res_data_g& resData;
};



void*
granularityCalculatingFunction(std::atomic<PageRankType>* pr_curr, std::atomic<PageRankType>* pr_next, void* data,
                               int thread_id) {
    timer serial_timer;
    serial_timer.start();

    arg_struct_g& args = *(arg_struct_g*) data;;
    Graph& g = args.graph;
    res_data_g& resData = args.resData;
    uintV max_iterations = resData.max_iterations;

    int granularity = resData.granularity;

    CustomBarrier& b = args.b;

    timer t;
    t.start();
    double barrier1_time = 0;
    double barrier2_time = 0;

    resData.edge_arr[thread_id] = 0;
    resData.vertices_processed[thread_id] = 0; // initialized to 0, c++ does not have default 0 value

    for (int i = 0; i < max_iterations; i++) {
        while (true) {
            uintE u = resData.getNextVertexToBeProcessed();
            if (u <= -1) break;
            for (uintV j = u; j < u + granularity; j++) {

                resData.vertices_processed[thread_id]++;

                uintE out_degree = g.vertices_[u].getOutDegree();
                for (uintE j = 0; j < out_degree; j++) {
                    uintV v = g.vertices_[u].getOutNeighbor(j);

                    float temp = pr_next[v];
                    while (!pr_next[v].compare_exchange_weak(temp, pr_next[v] + pr_curr[u] / out_degree));
                }
                if (u >= g.n_) break; // n is the total number of vertices in the graph
            }
        }
        b.wait();
        barrier1_time += t.stop();
        if (thread_id == 0)
            resData.starting_index = g.n_;

        b.wait();
        t.start();

        while (true) {
            uintE v = resData.getNextVertexToBeProcessed();
            if (v <= -1) break;

            for (uintV j = v; j < v + granularity; j++) {
                pr_next[v] = PAGE_RANK(pr_next[v]);
                // reset pr_curr for the next iteration
                pr_curr[v] = pr_next[v].load();
                pr_next[v] = 0.0;

                if(v >= g.n_) break; // n is the total number of vertices in the graph
            }
        }
        b.wait();
        barrier2_time += t.stop();
    }
}

void pageRankParallel(Graph& g, uint max_iterations, int n_workers, uint granularity) {
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

    res_data_g resData{};
    resData.max_iterations = max_iterations;
    resData.vertices_processed = new uintV[n_workers];
    resData.time_taken_s = new double[n_workers];
    resData.edge_arr = new int32_t[n_workers];
    resData.barrier1_time_arr = new int32_t[n_workers];
    resData.barrier2_time_arr = new int32_t[n_workers];

    resData.starting_index = g.n_ - 1;
    resData.granularity = granularity;

    arg_struct_g aStruct{resData, b, g};

    for (int i = 0; i < n_workers; i++)
        threads[i] = std::thread{granularityCalculatingFunction, pr_curr, pr_next, &aStruct, i};

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
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    delete[] pr_curr;
    delete[] pr_next;

    delete[] resData.vertices_processed;
    delete[] resData.time_taken_s;
    delete[] resData.edge_arr;
    delete[] resData.barrier1_time_arr;
    delete[] resData.barrier2_time_arr;
}

int main(int argc, char* argv[]) {
    cxxopts::Options options(
            "page_rank_push",
            "Calculate page_rank using serial and parallel execution");
    options.add_options(
            "",
            {
                    {"nWorkers",    "Number of workers",
                            cxxopts::value<uint>()->default_value(
                                    DEFAULT_NUMBER_OF_WORKERS)},
                    {"nIterations", "Maximum number of iterations",
                            cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                    {"granularity", "Granularity to be used", cxxopts::value<uint>()->default_value(
                            DEFAULT_GRANULARITY)},
                    {"inputFile",   "Input graph file path",
                            cxxopts::value<std::string>()->default_value(
                                    "/scratch/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    uint granularity = cl_options["granularity"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "Granularity : " << granularity << "\n";
    std::cout << "Iterations : " << max_iterations << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";


    switch (granularity) {
        case 0:
            break;
        case 1:
            pageRankParallel(g, max_iterations, n_workers);
            break;
        default:
            pageRankParallel(g, max_iterations, n_workers, granularity);
            break;
    }
    return 0;
}


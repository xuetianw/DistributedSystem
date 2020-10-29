//
// Created by fred on 2020-10-16.
//

#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <atomic>

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

std::mutex mtx;           // mutex for critical section


//struct barrier_type {
//    barrier_type(uint numOfWorker) : num_of_worker(numOfWorker) {}
//    int active_worker_count;
//    int num_of_worker;
//    std::condition_variable cv;
//
//    std::mutex lock;
//
//    void waits()
//    {
//        std::unique_lock<std::mutex> lk(lock);
////        std::cerr << "active_worker_count : " << active_worker_count << std::endl;
//
//        if (active_worker_count > 1) {
//            active_worker_count--;
//            std::cerr << "Waiting... \n";
//            cv.wait(lk);
//            std::cerr << "...finished waiting.\n";
//        }  else if (active_worker_count == 1){
//            cv.notify_all();
//            active_worker_count = num_of_worker;
//            std::cerr << "active_worker_count reset to : " << active_worker_count << std::endl;
//        }
//    }
//};

struct res_data {
    uintV* num_vertices_arr;
    double* time_taken_s;
    uintV max_iterations;
};

struct arg_struct {
    arg_struct(res_data& resData, CustomBarrier& b, Graph& graph) : resData(resData), graph(graph), b(b) {}

    Graph& graph;
//    barrier_type& b;dd
    CustomBarrier& b;
    res_data& resData;
};

std::mutex mutex;

void* calculatingFunction(std::atomic<PageRankType>* pr_curr, std::atomic<PageRankType>* pr_next, void* data, int thread_id) {

    timer serial_timer;
    serial_timer.start();

    arg_struct& args = *(arg_struct*) data;;
    Graph& g = args.graph;
    res_data& resData = args.resData;
    uintV max_iterations = resData.max_iterations;

    uint num_of_vertices = resData.num_vertices_arr[thread_id];

    uintV starting_vertex = thread_id == 0 ? 0 :
                            thread_id * resData.num_vertices_arr[thread_id - 1];

    uintV ending_vertex= starting_vertex + num_of_vertices - 1;

//    std::cerr << "starting_vertex: " << starting_vertex << "  ending_vertex : " << ending_vertex << std::endl;

//    barrier_type& b = args.b;
    CustomBarrier& b = args.b;


    for (int i = 0; i < max_iterations; i++) {
        for (uintV u = starting_vertex; u <= ending_vertex; u++) {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE j = 0; j < out_degree; j++) {
                uintV v = g.vertices_[u].getOutNeighbor(j);

                float temp = pr_next[v];
                while(!pr_next[v].compare_exchange_weak(temp, pr_next[v] + pr_curr[u] / out_degree))
                    ;
            }
        }
        b.wait();

        for (uintV v = starting_vertex; v <= ending_vertex; v++) {
            pr_next[v] = PAGE_RANK(pr_next[v]);
            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v].load();
            pr_next[v] = 0.0;
        }

        b.wait();
    }

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
    resData.num_vertices_arr = new uintV[n_workers];
    for (int i = 0; i < n_workers - 1; i++)
        resData.num_vertices_arr[i] = n / n_workers;

    resData.num_vertices_arr[n_workers - 1] =
            n % n_workers == 0 ? n / n_workers : n - (n / n_workers) * (n_workers - 1);
    resData.max_iterations = max_iterations;
    resData.time_taken_s = new double[n_workers];

    arg_struct aStruct{resData, b, g};

    for (int i = 0; i < n_workers; i++)
        threads[i] = std::thread{calculatingFunction, pr_curr, pr_next, &aStruct, i};

    std::cout << "thread_id, time_taken\n";

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
        printf("%d,\t %f\n", i, resData.time_taken_s[i]);
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

    delete[] resData.time_taken_s;
    delete[] resData.num_vertices_arr;
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
                    {"inputFile",   "Input graph file path",
                            cxxopts::value<std::string>()->default_value(
                                    "/scratch/input_graphs/roadNet-CA")},
            });

    auto cl_options = options.parse(argc, argv);
    uint n_workers = cl_options["nWorkers"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
//    uint n_workers = 4;
//    uint max_iterations = 10;
//    std::string input_file_path = "/home/fred/DistributedSystemProject/assignment1/inputfiles/lj";

#ifdef USE_INT
    std::cout << "Using INT\n";
#else
    std::cout << "Using FLOAT\n";
#endif
    std::cout << std::fixed;
    std::cout << "Number of workers : " << n_workers << "\n";
    std::cout << "max_iterations : " << max_iterations << "\n";



    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    pageRankVertexBasedParallel(g, max_iterations, n_workers);

//    std::cout << "Number of workers : \n" << n_workers << "\n";

    return 0;
}

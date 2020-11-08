#include "core/circular_queue.h"
#include "core/utils.h"
#include <pthread.h>
#include <stdlib.h>



struct res_data_struct {
    long n_items;
    int& active_producer_count;
    int& active_consumer_count;
    int puc_csm_id;

    long items_consumed;
    long items_produced;
    long value_produced;
    long value_consumed;
    double time_taken;
    CircularQueue& circularQueue;

    // Declaration of thread condition variable
    explicit res_data_struct(long n_items, CircularQueue& circularQueue, int id,
                             int& active_producer_count, int& active_consumer_count);

};

class ProducerConsumerProblem {
    long n_items;
    int n_producers;
    int n_consumers;
    CircularQueue production_buffer;

    pthread_t* producer_threads;
    pthread_t* consumer_threads;

    int active_producer_count;
    int active_consumer_count;



public:
    // The following 6 methods should be defined in the implementation file (solution.cpp)
    ProducerConsumerProblem(long _n_items, int _n_producers, int _n_consumers,
                            long _queue_size);

    ~ProducerConsumerProblem();

    void startProducers();

    void startConsumers();

    void joinProducers();

    void joinConsumers();

    void printStats();

    std::vector<res_data_struct> consumer_res_s;
    std::vector<res_data_struct> producer_res_s;
};

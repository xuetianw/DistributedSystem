#include "solution.h"

// Declaration of thread condition variable
pthread_cond_t cond1 = PTHREAD_COND_INITIALIZER;
pthread_cond_t cond2 = PTHREAD_COND_INITIALIZER;

pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t producer_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t consumer_mutex = PTHREAD_MUTEX_INITIALIZER;

void* producerFunction(void* _arg) {
    res_data_struct& res_data = *(res_data_struct*) _arg;
    long threshold = res_data.n_items;
    long data = 0;
    long items_produced = 0;
    long value_produced = 0;
    timer t;
    t.start();

    while (items_produced < threshold) {
        pthread_mutex_lock(&queue_mutex);
        {
            bool ret = res_data.circularQueue.enqueue(data);
            if (ret) {
                if (res_data.circularQueue.itemCount() == 1) {
                    pthread_cond_broadcast(&cond1);
                }
                value_produced += data;
                items_produced++;
                data++;
            } else {
                pthread_cond_wait(&cond2, &queue_mutex);
            }
        }
        pthread_mutex_unlock(&queue_mutex);
    }
    double time_taken = t.stop();
    res_data.items_produced = items_produced;
    res_data.value_produced = value_produced;
    res_data.time_taken = time_taken;

    pthread_mutex_lock(&producer_mutex);
    {
        res_data.active_producer_count--;
    }
    pthread_mutex_unlock(&producer_mutex);

    if (res_data.active_producer_count == 0)
        pthread_cond_broadcast(&cond1);

}

void* consumerFunction(void* _arg) {
    res_data_struct& res_data = *(res_data_struct*) _arg;
    CircularQueue& productionQueue = res_data.circularQueue;
    timer t;
    t.start();

    long item;
    long items_consumed = 0;
    long value_consumed = 0;
    while (true) {
        pthread_mutex_lock(&queue_mutex);
        {
            bool ret = productionQueue.dequeue(&item);
            if (ret) {
                if (productionQueue.itemCount() == productionQueue.getCapacity() - 1) {
                    pthread_cond_broadcast(&cond2);
                }
                value_consumed += item;
                items_consumed++;
            } else {
                if (res_data.active_producer_count == 0) {
                    pthread_mutex_unlock(&queue_mutex);
                    break;
                }
                pthread_cond_wait(&cond1, &queue_mutex);
            }
        }
        pthread_mutex_unlock(&queue_mutex);
    }

    double time_taken = t.stop();
    res_data.items_consumed = items_consumed;
    res_data.value_consumed = value_consumed;
    res_data.time_taken = time_taken;

    pthread_mutex_lock(&consumer_mutex);
    {
        res_data.active_consumer_count--;
    }
    pthread_mutex_unlock(&consumer_mutex);
}

ProducerConsumerProblem::ProducerConsumerProblem(long _n_items,
                                                 int _n_producers,
                                                 int _n_consumers,
                                                 long _queue_size)
        : n_items(_n_items), n_producers(_n_producers), n_consumers(_n_consumers), production_buffer(_queue_size) {

    active_producer_count = _n_producers;
    active_consumer_count = _n_consumers;

    for (int i = 0; i < n_producers; ++i)
        producer_res_s.emplace_back(_n_items, production_buffer, i, active_producer_count,
                                    active_consumer_count);

    for (int i = 0; i < n_consumers; ++i)
        consumer_res_s.emplace_back(_n_items, production_buffer, i, active_producer_count,
                                    active_consumer_count);

    std::cout << "Constructor\n";
    std::cout << "Number of items: " << n_items << "\n";
    std::cout << "Number of producers: " << n_producers << "\n";
    std::cout << "Number of consumers: " << n_consumers << "\n";
    std::cout << "Queue size: " << _queue_size << "\n";
    producer_threads = new pthread_t[n_producers];
    consumer_threads = new pthread_t[n_consumers];
}

ProducerConsumerProblem::~ProducerConsumerProblem() {
    std::cout << "Destructor\n";
    delete[] producer_threads;
    delete[] consumer_threads;
}


void ProducerConsumerProblem::startProducers() {
    for (int i = 0; i < n_producers; i++)
        pthread_create(&producer_threads[i], nullptr, producerFunction, &producer_res_s[i]);
}

void ProducerConsumerProblem::startConsumers() {
    for (int i = 0; i < n_consumers; i++)
        pthread_create(&consumer_threads[i], nullptr, consumerFunction, &consumer_res_s[i]);
}

void ProducerConsumerProblem::joinProducers() {
    // Join the producer threads with the main thread using pthread_join
    for (int i = 0; i < n_producers; i++)
        pthread_join(producer_threads[i], nullptr);
}

void ProducerConsumerProblem::joinConsumers() {
    // Join the consumer threads with the main thread using pthread_join
    for (int i = 0; i < n_consumers; i++)
        pthread_join(consumer_threads[i], nullptr);
}

void ProducerConsumerProblem::printStats() {
    std::cout << "Producer stats\n";
    std::cout << "producer_id, items_produced, value_produced, time_taken \n";

    long total_produced = 0;
    long total_value_produced = 0;
    for (int i = 0; i < n_producers; i++) {
        std::cout << i << ",\t" << producer_res_s[i].items_produced << ",\t"
                  << producer_res_s[i].value_produced << ",\t" << producer_res_s[i].time_taken << std::endl;

        total_produced += producer_res_s[i].items_produced;
        total_value_produced += producer_res_s[i].value_produced;
    }
    std::cout << "Total produced = " << total_produced << "\n";
    std::cout << "Total value produced = " << total_value_produced << "\n";
    std::cout << "Consumer stats\n";
    std::cout << "consumer_id, items_consumed, value_consumed, time_taken \n";

    long total_consumed = 0;
    long total_value_consumed = 0;
    for (int i = 0; i < n_consumers; i++) {
        res_data_struct csm = consumer_res_s[i];
        std::cout << i << ",\t" << csm.items_consumed << ",\t"
                  << csm.value_consumed << ",\t" << csm.time_taken << std::endl;

        total_consumed += csm.items_consumed;
        total_value_consumed += csm.value_consumed;
    }
    std::cout << "Total consumed = " << total_consumed << "\n";
    std::cout << "Total value consumed = " << total_value_consumed << "\n";
}


res_data_struct::res_data_struct(long n_items, CircularQueue& circularQueue, int id, int& active_producer_count,
                                 int& active_consumer_count) :
        circularQueue(circularQueue), puc_csm_id(id),
        active_producer_count(active_producer_count),
        active_consumer_count(active_consumer_count),
        n_items(n_items) {
    time_taken = 0;
    items_consumed = 0;
    value_produced = 0;
    items_produced = 0;
    value_consumed = 0;
}

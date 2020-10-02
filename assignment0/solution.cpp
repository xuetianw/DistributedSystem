#include "solution.h"

// Declaration of thread condition variable
pthread_cond_t cond1 = PTHREAD_COND_INITIALIZER;
pthread_cond_t cond2 = PTHREAD_COND_INITIALIZER;

pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t producer_mutex2 = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t consumer_mutex = PTHREAD_MUTEX_INITIALIZER;

void* producerFunction(void* _arg) {
    // Parse the _arg passed to the function.
    // Enqueue `n` items into the `production_buffer`. The items produced should
    // be 0, 1, 2,..., (n-1).
    // Keep track of the number of items produced and the value produced by the
    // thread.
    // The producer that was last active should ensure that all the consumers have
    // finished. NOTE: Each thread will enqueue `n` items.
    // Use mutex variables and conditional variables as necessary.
    // Each producer enqueues `n` items

    res_data_struct& res_data = *(res_data_struct*) _arg;

    long threshold = res_data.n_items;
    long item = 0;
    long items_produced = 0;
    long value_produced = 0;
    timer t;
    t.start();


    while (items_produced < threshold) {

        pthread_mutex_lock(&queue_mutex);

        bool ret = res_data.circularQueue.enqueue(item);

        if (ret) {

            if (res_data.circularQueue.itemCount() == 1) {
                // The queue is no longer empty
                // Signal all consumers indicating queue is not empty
                pthread_cond_broadcast(&cond2);
                printf("signal cond2 \n");
            }
//            pthread_mutex_lock(&mutex);

            value_produced += item;
            items_produced++;
            item++;

        } else {
            // production_buffer is full, so block on conditional variable waiting for consumer to signal.
            printf("buffer full \n");
            pthread_cond_wait(&cond1, &queue_mutex);
            printf("producerFunction awaken because buffer not empty \n");

//            res_data.circularQueue.enqueue(item);
//            value_produced += item;
//            items_produced++;
//            item++;

        }
        pthread_mutex_unlock(&queue_mutex);
    }
    double time_taken = t.stop();
    res_data.items_produced = items_produced;
    res_data.value_produced = value_produced;
    res_data.time_taken = time_taken;
    // After production is completed:
    // Update the number of producers that are currently active.

    pthread_mutex_lock(&producer_mutex2);

    res_data.active_producer_count--;

    pthread_mutex_unlock(&producer_mutex2);
//    printf("producer work done \n");
    // The producer that was last active (can be determined using `active_producer_count`) will keep signalling the consumers until all consumers have finished (can be determined using `active_consumer_count`).
    if (res_data.active_producer_count == 0) {
        printf("last producer work done \n");
        printf("items_produced value is now: %ld\n", items_produced);
        pthread_cond_broadcast(&cond2);
    }
    printf("producer done\n");
}

void* consumerFunction(void* _arg) {
    // Parse the _arg passed to the function.
    // The consumer thread will consume items by dequeueing the items from the
    // `production_buffer`.
    // Keep track of the number of items consumed and the value consumed by the
    // thread.
    // Once the productions is complete and the queue is also empty, the thread
    // will exit. NOTE: The number of items consumed by each thread need not be
    // same.
    // Use mutex variables and conditional variables as necessary.

    // Each consumer dequeues items from the `production_buffer`


    res_data_struct& res_data = *(res_data_struct*) _arg;

    CircularQueue& productionQueue = res_data.circularQueue;
    timer t;
    t.start();

    long item;
    long items_consumed = 0;
    long value_consumed = 0;
    while (true) {

        pthread_mutex_lock(&queue_mutex);

        bool ret = productionQueue.dequeue(&item);

//        pthread_mutex_unlock(&queue_mutex);

        if (ret) {

            if (productionQueue.itemCount() == productionQueue.getCapacity() - 1) {
                // The queue is no longer full
                // Signal all producers indicating
                // is not full
                pthread_cond_broadcast(&cond1);
                printf("signal cond1 \n");
            }

//            pthread_mutex_lock(&mutex1);

            value_consumed += item;
            items_consumed++;

//            pthread_mutex_unlock(&mutex1);

        } else {
            // production_buffer is empty, so block on conditional variable waiting for producer to signal.
            // The thread can wake up because of 2 scenarios:
            // Scenario 1: There are no more active producers (i.e., production is complete) and the queue is empty.
            // This is the exit condition for consumers, and at this point consumers should decrement `active_consumer_count`.
            // Scenario 2: The queue is not empty and/or the producers are active. Continue consuming.

            printf("buffer empty \n");

//            if (!res_data.circularQueue.isEmpty() && res_data.active_producer_count == 0) {
//                while (productionQueue.dequeue(&item)) {
//                    value_consumed += item;
//                    items_consumed++;
//                }
//                break;
//            }


            if (res_data.active_producer_count != 0) {
                printf("11111111111111111 \n");
                pthread_cond_wait(&cond2, &queue_mutex);
            }


            printf("consumerFunction awaken because buffer is no longer empty or last buffer producer work done \n");

            if (res_data.active_producer_count == 0 && res_data.circularQueue.isEmpty()) {
                printf("test\n");
                pthread_mutex_unlock(&queue_mutex);
                break;
            }

//            if (!res_data.circularQueue.isEmpty() || res_data.active_producer_count != 0) {
//
//                productionQueue.dequeue(&item);
//                value_consumed += item;
//                items_consumed++;
//
//            }
//            printf("continue consuming\n");
        }
        pthread_mutex_unlock(&queue_mutex);
    }



    double time_taken = t.stop();
    res_data.items_consumed = items_consumed;
    res_data.value_consumed = value_consumed;
    res_data.time_taken = time_taken;


    pthread_mutex_lock(&consumer_mutex);

    res_data.active_consumer_count--;

    pthread_mutex_unlock(&consumer_mutex);

    printf("buffer consumer work done \n");
}

ProducerConsumerProblem::ProducerConsumerProblem(long _n_items,
                                                 int _n_producers,
                                                 int _n_consumers,
                                                 long _queue_size)
        : n_items(_n_items), n_producers(_n_producers), n_consumers(_n_consumers), production_buffer(_queue_size) {


    for (int i = 0; i < n_producers; ++i) {
        producer_res_s.push_back(
                new res_data_struct(_n_items, production_buffer, i, active_producer_count,
                                    active_consumer_count,cond1));
    }

    for (int i = 0; i < n_consumers; ++i) {
        consumer_res_s.push_back(
                new res_data_struct(_n_items, production_buffer, i, active_producer_count,
                                    active_consumer_count,cond1));
    }

    std::cout << "Constructor\n";
    std::cout << "Number of items: " << n_items << "\n";
    std::cout << "Number of producers: " << n_producers << "\n";
    std::cout << "Number of consumers: " << n_consumers << "\n";
    std::cout << "Queue size: " << _queue_size << "\n";
    producer_threads = new pthread_t[n_producers];
    consumer_threads = new pthread_t[n_consumers];

    active_producer_count = _n_producers;
    active_consumer_count = _n_consumers;
    printf("active_producer_count %d \n", active_producer_count);
    printf("active_consumer_count %d \n", active_consumer_count);
    printf("\n");
    // Initialize all mutex and conditional variables here.
}

ProducerConsumerProblem::~ProducerConsumerProblem() {
    std::cout << "Destructor\n";
    delete[] producer_threads;
    delete[] consumer_threads;
    // Destroy all mutex and conditional variables here.
}


void* PrintHello(void* threadid) {
    long tid;
    tid = (long) threadid;
    std::cout << "Hello World! Thread ID, " << tid << std::endl;
    pthread_exit(NULL);
}


void ProducerConsumerProblem::startProducers() {
//    std::cout << "Starting Producers\n";
    active_producer_count = n_producers;
    // Create producer threads P1, P2, P3,.. using pthread_create.
    for (int i = 0; i < n_producers; i++) {
        pthread_create(&consumer_threads[i], nullptr, producerFunction, producer_res_s[i]);
    }
}

void ProducerConsumerProblem::startConsumers() {
    std::cout << "Starting Consumers\n";
    active_consumer_count = n_consumers;
    // Create consumer threads C1, C2, C3,.. using pthread_create.

    for (int i = 0; i < n_consumers; i++) {
        pthread_create(&producer_threads[i], nullptr, consumerFunction, consumer_res_s[i]);
        printf("consumer number : %d\n", i);
    }
}

void ProducerConsumerProblem::joinProducers() {
    std::cout << "Joining Producers\n";
    // Join the producer threads with the main thread using pthread_join
    for (int i = 0; i < n_producers; i++) {
        pthread_join(producer_threads[i], (void**) producer_res_s[i]);
        printf("done Joining Producers %d \n", i);
    }
    std::cout << "done Joining Producers\n";
}

void ProducerConsumerProblem::joinConsumers() {
    std::cout << "Joining Consumers\n";
    // Join the consumer threads with the main thread using pthread_join
    for (int i = 0; i < n_consumers; i++) {
        pthread_join(consumer_threads[i], (void**) consumer_res_s[i]);
        printf("done Joining Consumers %d \n", i);
    }

    std::cout << "done Joining Consumers\n";
}

void ProducerConsumerProblem::printStats() {
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "Producer stats\n";
    std::cout << "producer_id, items_produced, value_produced, time_taken \n";
    // Make sure you print the producer stats in the following manner
    // 0, 1000, 499500, 0.00123596
    // 1, 1000, 499500, 0.00154686
    // 2, 1000, 499500, 0.00122881
    long total_produced = 0;
    long total_value_produced = 0;
    for (int i = 0; i < n_producers; i++) {
        // ---
        //
        // ---
        std::cout << i << "\t\t\t" << producer_res_s[i]->items_produced << "\t\t\t"
                << producer_res_s[i]->value_produced << "\t\t\t" << producer_res_s[i]->time_taken << std::endl;

        total_produced += producer_res_s[i]->items_produced;
        total_value_produced += producer_res_s[i]->value_produced;
    }
    std::cout << "Total produced = " << total_produced << "\n";
    std::cout << "Total value produced = " << total_value_produced << "\n";
    std::cout << "Consumer stats\n";
    std::cout << "consumer_id, items_consumed, value_consumed, time_taken \n";
    // Make sure you print the consumer stats in the following manner
    // 0, 677, 302674, 0.00147414
    // 1, 648, 323301, 0.00142694
    // 2, 866, 493382, 0.00139689
    // 3, 809, 379143, 0.00134516
    long total_consumed = 0;
    long total_value_consumed = 0;
    for (int i = 0; i < n_consumers; i++) {
        // ---
        //
        // ---
        std::cout << i << "\t\t\t" << consumer_res_s[i]->items_consumed << "\t\t\t"
                << consumer_res_s[i]->value_consumed << "\t\t\t" << consumer_res_s[i]->time_taken << std::endl;

        total_consumed += consumer_res_s[i]->items_consumed;
        total_value_consumed += consumer_res_s[i]->value_consumed;
    }
    std::cout << "Total consumed = " << total_consumed << "\n";
    std::cout << "Total value consumed = " << total_value_consumed << "\n";
}



res_data_struct::res_data_struct(long n_items, CircularQueue& circularQueue, int id, int& active_producer_count,
                                 int& active_consumer_count, pthread_cond_t& pthread_cond_t) :
        circularQueue(circularQueue), puc_csm_id(id),
        active_producer_count(active_producer_count),
        active_consumer_count(active_consumer_count),
        n_items(n_items), cond1(pthread_cond_t) {
    time_taken = 0;
    items_consumed = 0;
    value_produced = 0;
}

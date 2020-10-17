//
// Created by fred on 2020-10-15.
//

//#include <atomic>
//template<typename T>
//struct node
//{
//    T data;
//    node* next;
//    node(const T& data) : data(data), next(nullptr) {}
//};
//
//template<typename T>
//class stack
//{
//    std::atomic<node<T>*> head;
//public:
//    void push(const T& data)
//    {
//        node<T>* new_node = new node<T>(data);
//
//        // put the current value of head into new_node->next
//        new_node->next = head.load(std::memory_order_relaxed);
//
//        // now make new_node the new head, but if the head
//        // is no longer what's stored in new_node->next
//        // (some other thread must have inserted a node just now)
//        // then put that new head into new_node->next and try again
//        while(!head.compare_exchange_weak(new_node->next, new_node,
//                                          std::memory_order_release,
//                                          std::memory_order_relaxed))
//            ; // the body of the loop is empty
//
//// Note: the above use is not thread-safe in at least
//// GCC prior to 4.8.3 (bug 60272), clang prior to 2014-05-05 (bug 18899)
//// MSVC prior to 2014-03-17 (bug 819819). The following is a workaround:
////      node<T>* old_head = head.load(std::memory_order_relaxed);
////      do {
////          new_node->next = old_head;
////       } while(!head.compare_exchange_weak(old_head, new_node,
////                                           std::memory_order_release,
////                                           std::memory_order_relaxed));
//    }
//};
//int main()
//{
//    stack<int> s;
//    s.push(1);
//    s.push(2);
//    s.push(3);
//}



//#include <string>
//#include <thread>
//#include <vector>
//#include <iostream>
//#include <atomic>
//#include <chrono>
//
//// meaning of cnt:
//// 5: there are no active readers or writers.
//// 1...4: there are 4...1 readers active, The writer is blocked
//// 0: temporary value between fetch_sub and fetch_add in reader lock
//// -1: there is a writer active. The readers are blocked.
//const int N = 5; // four concurrent readers are allowed
//std::atomic<int> cnt(N);
//
//std::vector<int> data;
//
//void reader(int id)
//{
//    for(;;)
//    {
//        // lock
//        while(std::atomic_fetch_sub(&cnt, 1) <= 0)
//            std::atomic_fetch_add(&cnt, 1);
//        // read
//        if(!data.empty())
//            std::cout << (  "reader " + std::to_string(id)
//                            + " sees " + std::to_string(*data.rbegin()) + '\n');
//        if(data.size() == 25)
//            break;
//        // unlock
//        std::atomic_fetch_add(&cnt, 1);
//        // pause
//        std::this_thread::sleep_for(std::chrono::milliseconds(1));
//    }
//}
//
//void writer()
//{
//    for(int n = 0; n < 25; ++n)
//    {
//        // lock
//        while(std::atomic_fetch_sub(&cnt, N+1) != N)
//            std::atomic_fetch_add(&cnt, N+1);
//        // write
//        data.push_back(n);
//        std::cout << "writer pushed back " << n << '\n';
//        // unlock
//        std::atomic_fetch_add(&cnt, N+1);
//        // pause
//        std::this_thread::sleep_for(std::chrono::milliseconds(1));
//    }
//}
//
//int main()
//{
//    std::vector<std::thread> v;
//    for (int n = 0; n < N; ++n) {
//        v.emplace_back(reader, n);
//    }
//    v.emplace_back(writer);
//    for (auto& t : v) {
//        t.join();
//    }
//}


//#include <atomic>
//template<typename T>
//struct node
//{
//    T data;
//    node* next;
//    node(const T& data) : data(data), next(nullptr) {}
//};
//
//template<typename T>
//class stack
//{
//    std::atomic<node<T>*> head;
//public:
//    void push(const T& data)
//    {
//        node<T>* new_node = new node<T>(data);
//
//        // put the current value of head into new_node->next
//        new_node->next = head.load(std::memory_order_relaxed);
//
//        // now make new_node the new head, but if the head
//        // is no longer what's stored in new_node->next
//        // (some other thread must have inserted a node just now)
//        // then put that new head into new_node->next and try again
//        while(!head.compare_exchange_weak(new_node->next, new_node,
//                                          std::memory_order_release,
//                                          std::memory_order_relaxed))
//            ; // the body of the loop is empty
//
//// Note: the above use is not thread-safe in at least
//// GCC prior to 4.8.3 (bug 60272), clang prior to 2014-05-05 (bug 18899)
//// MSVC prior to 2014-03-17 (bug 819819). The following is a workaround:
////      node<T>* old_head = head.load(std::memory_order_relaxed);
////      do {
////          new_node->next = old_head;
////       } while(!head.compare_exchange_weak(old_head, new_node,
////                                           std::memory_order_release,
////                                           std::memory_order_relaxed));
//    }
//};
//int main()
//{
//    stack<int> s;
//    s.push(1);
//    s.push(2);
//    s.push(3);
//}


#include <atomic>


//template<typename T>
class stack {
    std::atomic<float> head;

public:
    void add(const float data)
    {
//        node<T>* new_node = new node<T>(data);

        // put the current value of head into new_node->next
//        new_node->next = head.load(std::memory_order_relaxed);

        // now make new_node the new head, but if the head
        // is no longer what's stored in new_node->next
        // (some other thread must have inserted a node just now)
        // then put that new head into new_node->next and try again
        float temp = head;
        while(!head.compare_exchange_weak(temp, temp + data,
                                          std::memory_order_release,
                                          std::memory_order_relaxed))
            ; // the body of the loop is empty

// Note: the above use is not thread-safe in at least
// GCC prior to 4.8.3 (bug 60272), clang prior to 2014-05-05 (bug 18899)
// MSVC prior to 2014-03-17 (bug 819819). The following is a workaround:
//      node<T>* old_head = head.load(std::memory_order_relaxed);
//      do {
//          new_node->next = old_head;
//       } while(!head.compare_exchange_weak(old_head, new_node,
//                                           std::memory_order_release,
//                                           std::memory_order_relaxed));
    }
};
int main()
{
//    stack<int> s;
//    s.push(1);
//    s.push(2);
//    s.push(3);
}
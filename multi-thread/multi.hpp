#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>

class ThreadPool {
public:
    // Constructor to initialize the thread pool with a given number of threads
    explicit ThreadPool(size_t numThreads) :Thread_num(numThreads), stop(false), busy(0) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        // Acquire lock to access the task queue
                        std::unique_lock<std::mutex> lock(queueMutex);
                        // Wait until there are tasks or the pool is stopping
                        condition.wait(lock, [this] { return stop || !tasks.empty(); });
                        if (stop && tasks.empty()) return; // Exit thread if stopping and no tasks
                        ++busy;
                        // Get the next task from the queue
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    // Execute the task
                    task();
                    {
                        // Lock to update the busy count and notify finished condition
                        std::lock_guard<std::mutex> lock(queueMutex);
                        --busy;
                        cv_finished.notify_one();
                    }
                }
            });
        }
    }

    // Template function to execute a given task with arguments
    template<typename Func, typename... Args>
    void execute(Func&& func, Args&&... args) {
        // Wait until there's enough capacity in the queue
        while (tasks.size() > workers.size())
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        {
            // Lock and add the task to the queue
            std::lock_guard<std::mutex> lock(queueMutex);
            tasks.emplace(std::bind(std::forward<Func>(func), std::forward<Args>(args)...));
        }
        // Notify one worker thread to pick up the new task
        condition.notify_one();
    }

    // Overloaded execute function for member function pointers
    template<typename Func, typename... Args>
    void execute(void (Func::* func)(Args...), Func* obj, Args&&... args) {
        // Wait until there's enough capacity in the queue
        while (tasks.size() > workers.size())
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        {
            // Lock and add the member function task to the queue
            std::lock_guard<std::mutex> lock(queueMutex);
            tasks.emplace(std::bind(std::mem_fn(func), obj, std::forward<Args>(args)...));
        }
        // Notify one worker thread to pick up the new task
        condition.notify_one();
    }

    // Wait until all tasks in the queue are completed
    void waitFinished() {
        std::unique_lock<std::mutex> lock(queueMutex);
        cv_finished.wait(lock, [this] { return tasks.empty() && (busy == 0); });
    }

    // Destructor to stop all threads and clean up
    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            stop = true; // Set the stop flag to true
        }
        condition.notify_all(); // Notify all worker threads to finish
        for (std::thread& worker : workers) {
            if (worker.joinable()) worker.join(); // Join all worker threads
        }
    }

    size_t Thread_num; // Number of threads in the pool
    std::vector<std::thread> workers; // Vector of worker threads
    std::mutex mutex_ns, mutex_fp; // Mutex for potential future use (currently unused)

private:
    std::queue<std::function<void()>> tasks; // Task queue
    std::mutex queueMutex; // Mutex for task queue synchronization
    std::condition_variable condition; // Condition variable for task availability
    std::condition_variable cv_finished; // Condition variable to signal when tasks are finished
    bool stop; // Flag to stop all threads
    size_t busy; // Count of currently busy threads
};

// External pointer to a ThreadPool instance
extern ThreadPool* threadPool0;

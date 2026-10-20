#include "gcrank/detail/thread_pool.h"

namespace gcrank {

ThreadPool::ThreadPool(unsigned num_threads) {
  if (num_threads < 1) num_threads = 1;
  workers_.reserve(num_threads);
  for (unsigned i = 0; i < num_threads; ++i)
    workers_.emplace_back([this] { worker_loop(); });
}

ThreadPool::~ThreadPool() {
  try {
    stop();
  } catch (...) {
    // A destructor must not throw.
  }
}

void ThreadPool::submit(std::move_only_function<void()> task) {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (stopping_ || abort_) return;
    tasks_.push_back(std::move(task));
  }
  cv_task_.notify_one();
}

void ThreadPool::worker_loop() {
  for (;;) {
    std::move_only_function<void()> task;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_task_.wait(lock, [this] { return !tasks_.empty() || stopping_ || abort_; });
      if (abort_ || tasks_.empty()) {
        // stopping_ with an empty queue, or an error: leave.
        return;
      }
      task = std::move(tasks_.front());
      tasks_.pop_front();
      ++active_;
    }
    try {
      task();
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!error_) error_ = std::current_exception();
      abort_ = true;
      tasks_.clear();
      --active_;
      cv_task_.notify_all();
      cv_done_.notify_all();
      return;
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      --active_;
      if (tasks_.empty() && active_ == 0) cv_done_.notify_all();
    }
  }
}

void ThreadPool::wait() {
  std::exception_ptr error;
  {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_done_.wait(lock, [this] { return abort_ || (tasks_.empty() && active_ == 0); });
    error = error_;
  }
  if (error) std::rethrow_exception(error);
}

void ThreadPool::stop() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    stopping_ = true;
  }
  cv_task_.notify_all();
  for (std::thread& worker : workers_)
    if (worker.joinable()) worker.join();
  workers_.clear();
  std::exception_ptr error;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    error = error_;
  }
  if (error) std::rethrow_exception(error);
}

}  // namespace gcrank

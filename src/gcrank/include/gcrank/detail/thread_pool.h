#pragma once
//
// Fixed-size worker pool with a std::deque task queue.
//
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace gcrank {

// A pool of `num_threads` workers (at least 1) draining a FIFO queue.
//
// Exceptions thrown by a task are captured in the slot that ran it and
// rethrown from wait()/stop() on the thread that owns the pool. Once a task
// throws, the pool aborts: queued tasks are dropped and submit() becomes a
// no-op. Long-running tasks are supported (a task may loop until the caller
// signals it to finish).
class ThreadPool {
public:
  explicit ThreadPool(unsigned num_threads);
  ~ThreadPool();

  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  unsigned size() const { return static_cast<unsigned>(workers_.size()); }

  // Enqueues a task; ignored once the pool is stopping or aborted.
  void submit(std::move_only_function<void()> task);

  // Blocks until all submitted tasks have finished (or the pool aborted) and
  // rethrows the first task exception, if any.
  void wait();

  // Signals the workers to exit after draining the queue, joins them, and
  // rethrows a task exception, if any. Called implicitly by the destructor.
  void stop();

private:
  void worker_loop();

  std::mutex mutex_;
  std::condition_variable cv_task_;
  std::condition_variable cv_done_;
  std::deque<std::move_only_function<void()>> tasks_;
  std::vector<std::thread> workers_;
  std::size_t active_ = 0;
  bool stopping_ = false;
  bool abort_ = false;
  std::exception_ptr error_;
};

}  // namespace gcrank

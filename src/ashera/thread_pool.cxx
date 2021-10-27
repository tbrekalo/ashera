#include "ashera/thread_pool.hpp"

namespace ashera {

namespace detail {

std::shared_ptr<thread_pool::ThreadPool> thread_pool;

}

auto InitThreadPool(std::uint32_t const n_threads) -> void {
  if (detail::thread_pool) {
    throw std::runtime_error(
        "[ashera::InitThreadPool] Invlaid attempt to reinitialize a global "
        "thread pool");
  }

  detail::thread_pool = std::make_shared<thread_pool::ThreadPool>(n_threads);
}

auto GetThreadPoolPtr() -> std::shared_ptr<thread_pool::ThreadPool> {
  if (!detail::thread_pool) {
    throw std::runtime_error(
        "[ashera::GetThreadPoolPtr] Global thread pool not initialized");
  }

  return detail::thread_pool;
}

}  // namespace ashera

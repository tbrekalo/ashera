#ifndef ASHERA_THREAD_POOL_
#define ASHERA_THREAD_POOL_

#include <cstdint>
#include <memory>
#include <stdexcept>

#include "thread_pool/thread_pool.hpp"

namespace ashera {

/**
 * @brief Initialize global thread pool instance (not thread safe init)
 */
auto InitThreadPool(std::uint32_t const n_threads) -> void;

/**
 * @brief Get a copy of std::shared_ptr<thread_pool::ThreadPool>
 */
auto GetThreadPoolPtr() -> std::shared_ptr<thread_pool::ThreadPool>;

}  // namespace ashera

#endif /* ASHERA_THREAD_POOL_ */

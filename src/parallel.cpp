#include "parallel.hpp"

#include <exception>
#include <thread>
#include <vector>

namespace kbsrc::parallel {

int run(int count, Task task, void* context) {
    std::atomic<bool> stop(false);
    std::vector<std::exception_ptr> errors(count > 1 ? count : 1);
    auto body = [&](int i) {
        try {
            task(context, stop);
        } catch (...) {
            errors[i] = std::current_exception();
            stop = true;
        }
    };
    std::vector<std::thread> pool;
    try {
        pool.reserve(errors.size() - 1);
        for (int i = 1; i < count; i++) pool.emplace_back(body, i);
    } catch (...) {
        // Fewer threads than asked for: the ones started and this one share the work.
    }
    body(0);
    for (auto& t : pool) t.join();
    for (auto& e : errors)
        if (e) std::rethrow_exception(e);
    return (int)pool.size() + 1;
}

}  // namespace kbsrc::parallel

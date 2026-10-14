// Worker threads for the engine copies. Compiled once, for the baseline instruction set (src/parallel.cpp).
#pragma once

#include <atomic>

namespace kbsrc::parallel {

using Task = void (*)(void* context, const std::atomic<bool>& stop);

// Runs task(context, stop) on up to `count` threads, the calling one included, and returns the number of threads
// that ran once all of them have finished. Each call must take its work items from a shared counter until they
// run out or `stop` is set, so that the result does not depend on how many threads take part; threads the system
// cannot start are left out. A call that throws (std::bad_alloc, say) sets `stop`; once every thread has
// finished, its exception is rethrown (the calling thread's first if several threw, then in start order).
int run(int count, Task task, void* context);

// Same for a callable f(stop).
template <class F> int run(int count, F& f) {
    return run(count, [](void* p, const std::atomic<bool>& stop) { (*static_cast<F*>(p))(stop); }, &f);
}

}  // namespace kbsrc::parallel

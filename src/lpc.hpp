// Burg LPC extrapolation used to extend the signal past its edges
#pragma once

#include <vector>

namespace kbsrc::detail {

// Predicts n samples following x[0..len) from an order-`order` Burg model trained on the last `train` samples.
// Returns zeros if the segment is silent, too short, or the prediction is unstable.
std::vector<double> burg_extrap(const double* x, long long len, long long n, int order = 32, int train = 4096);

}  // namespace kbsrc::detail

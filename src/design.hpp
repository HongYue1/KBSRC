// Kernel design: Kaiser-windowed sinc and its spectrum
#pragma once

#include <vector>

#include "dsp/double_double.hpp"

namespace kbsrc::detail {

using dsp::dd;

// One-sided taps h[0..L] of a Kaiser-windowed sinc with relative cutoff a/b (cutoff = (a/b) * fk/2),
// computed in long double with exact integer reduction of the sinc argument.
std::vector<dd> kernel_taps(int L, long long a, long long b, double beta);

// Same, with L = floor(T * fk).
std::vector<dd> kernel(long long fk, long long a, long long b, double T, double beta);

// Real spectrum (bins 0..n/2) of the symmetric kernel placed circularly in length n, divided by scale.
std::vector<dd> kernel_spectrum(const std::vector<dd>& h, int n, double scale);

}  // namespace kbsrc::detail

// Kernel design: Kaiser-windowed sinc and its spectrum
#pragma once

#include <vector>

#include "dsp/double_double.hpp"

namespace kbsrc::KBSRC_ISA::detail {

using dsp::dd;

// One-sided taps h[0..L] of a Kaiser-windowed sinc with relative cutoff a/b (cutoff = (a/b) * fk/2),
// computed in double-double with exact integer reduction of the sinc argument.
std::vector<dd> kernel_taps(int L, long long a, long long b, double beta);

// Real spectrum (bins 0..n/2) of the symmetric kernel placed circularly in length n, divided by scale.
std::vector<dd> kernel_spectrum(const std::vector<dd>& h, int n, double scale);

}  // namespace kbsrc::KBSRC_ISA::detail

// The engine (FFT, kernel design, block processing) is compiled once per instruction set and selected at
// run time. KBSRC_ISA names the namespace of each copy so the copies never share symbols.
#pragma once
#ifndef KBSRC_ISA
#define KBSRC_ISA base
#endif

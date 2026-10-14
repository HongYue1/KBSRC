"""No AVX instructions outside the AVX2 engine of a portable x86-64 build.

A portable build must start and run on any x86-64 CPU, so only the AVX2 engine, which is entered after a
run-time CPU check, may contain AVX (VEX- or EVEX-encoded) instructions. Everything else runs on every CPU:
start-up code (static initialisers run before main), the generic engine, the WAV writer and the library
functions that both engines share. The test disassembles the executable with objdump and fails for every such
function that contains an AVX instruction. Exits with 77 (skipped) without objdump, symbols or an AVX2 engine.
Standard library only.

Usage: test_isa.py --bin <kbsrc> [--objdump <objdump>]
"""
import argparse
import re
import shutil
import subprocess
import sys

SKIP = 77
FUNCTION = re.compile(r"^[0-9a-f]+ <(.+)>:$")
INSTRUCTION = re.compile(r"^\s*[0-9a-f]+:\s+(\S+)\s*(.*)$")
# Mnemonics that start with "v" without being AVX: VMX/SVM virtualisation and segment verification.
NOT_AVX = re.compile(r"^(verr|verw|vm(call|clear|func|launch|load|mcall|ptrld|ptrst|read|resume|run|save|write|xoff|xon))")
AVX2_ENGINE = "kbsrc::avx2::"
AVX2_MANGLED = "5kbsrc4avx2"  # kbsrc::avx2 in a mangled name objdump cannot demangle (e.g. C++20 constraints)
# Static initialisers run before main on every CPU, whichever translation unit they belong to.
STARTUP = re.compile(r"^(_GLOBAL__|__static_initialization|__cxx_global_var_init|global (con|de)structors keyed to)")


def in_avx2_engine(name):
    if STARTUP.match(name):
        return False
    return AVX2_ENGINE in name or (name.startswith("_Z") and AVX2_MANGLED in name)


def avx_functions(lines):
    """Map each function name to the number of AVX instructions in it."""
    found, name = {}, None
    for line in lines:
        m = FUNCTION.match(line)
        if m:
            name = m.group(1)
            continue
        m = INSTRUCTION.match(line)
        if not m or name is None:
            continue
        mnemonic, operands = m.groups()
        if (mnemonic.startswith("v") and not NOT_AVX.match(mnemonic)) or re.search(r"%[yz]mm|%k[0-7]\b", operands):
            found[name] = found.get(name, 0) + 1
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="path to the kbsrc executable")
    ap.add_argument("--objdump", default="", help="objdump to use (default: objdump on PATH)")
    args = ap.parse_args()
    objdump = args.objdump or shutil.which("objdump")
    if not objdump or not shutil.which(objdump):
        print("skipped: objdump not found (pass --objdump)")
        sys.exit(SKIP)
    r = subprocess.run([objdump, "-d", "-C", "--no-show-raw-insn", args.bin], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"objdump failed ({r.returncode}): {r.stderr.strip()}")
    lines = r.stdout.splitlines()
    names = [m.group(1) for m in map(FUNCTION.match, lines) if m]
    if not any(in_avx2_engine(n) for n in names):
        print(f"skipped: no {AVX2_ENGINE} functions in {args.bin} (not a portable x86-64 build, or stripped)")
        sys.exit(SKIP)
    bad = {n: c for n, c in avx_functions(lines).items() if not in_avx2_engine(n)}
    for n, c in sorted(bad.items()):
        print(f"FAIL  {c} AVX instruction(s) outside the AVX2 engine: {n[:200]}")
    print(f"{len(names)} functions checked, {len(bad)} with AVX instructions outside the AVX2 engine")
    if bad:
        raise SystemExit(f"{len(bad)} function(s) would fault on CPUs without AVX")
    print("all checks passed")


if __name__ == "__main__":
    main()

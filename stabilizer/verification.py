#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from pathlib import Path

try:
    import stim
except ImportError:
    print("Error: This script requires the 'stim' package. Install via: pip install stim", file=sys.stderr)
    sys.exit(1)


def build_random_stim_circuit(n_qubits: int, rounds: int, ops_per_round: int, seed: int):
    rng = stim.PRNG(seed)
    c = stim.Circuit()
    meas_order = []  # list of (qubit_index) in order of M ops
    for _ in range(rounds):
        for _ in range(ops_per_round):
            gate = rng.randint(5)  # 0:H,1:S,2:CX,3:RESET,4:M
            if gate == 2 and n_qubits >= 2:
                # CX with random distinct ctrl/target
                ctrl = rng.randint(n_qubits)
                target = rng.randint(n_qubits - 1)
                if target >= ctrl:
                    target += 1
                c.append_operation("CX", [ctrl, target])
            elif gate == 0:
                c.append_operation("H", [rng.randint(n_qubits)])
            elif gate == 1:
                c.append_operation("S", [rng.randint(n_qubits)])
            elif gate == 3:
                c.append_operation("R", [rng.randint(n_qubits)])  # reset to |0>
            else:
                q = rng.randint(n_qubits)
                c.append_operation("M", [q])
                meas_order.append(q)
    return c, len(meas_order)


def write_openqasm(qasm_path: Path, n_qubits: int, stim_circuit: stim.Circuit):
    # Convert the stim circuit ops to OpenQASM 2.0 text with measure ordering preserved.
    # We assume qelib1.inc supports h, s, cx, reset, measure.
    ops = []
    meas_count = 0
    for inst in stim_circuit:
        name = inst.name
        if name == "H":
            (q,) = map(int, inst.targets_copy())
            ops.append(f"h q[{q}];")
        elif name == "S":
            (q,) = map(int, inst.targets_copy())
            ops.append(f"s q[{q}];")
        elif name == "CX":
            q0, q1 = map(int, inst.targets_copy())
            ops.append(f"cx q[{q0}],q[{q1}];")
        elif name == "R":
            (q,) = map(int, inst.targets_copy())
            ops.append(f"reset q[{q}];")
        elif name == "M":
            (q,) = map(int, inst.targets_copy())
            ops.append(f"measure q[{q}] -> c[{meas_count}];")
            meas_count += 1
        else:
            raise ValueError(f"Unsupported Stim op in exporter: {name}")

    header = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"qreg q[{n_qubits}];",
        f"creg c[{meas_count}];",
    ]
    qasm_text = "\n".join(header + ops) + "\n"
    qasm_path.write_text(qasm_text)
    return meas_count


def run_cpp_backend(exe: Path, backend: str, qasm_path: Path, seed: int, meas_count: int):
    # Expects the binary to print a single line like: MEASURE: 0 1 1 0 ...
    cmd = [
        str(exe),
        "--backend", backend,
        "--qasm", str(qasm_path),
        "--seed", str(seed),
        "--meas-count", str(meas_count),
    ]
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[{backend}] runner failed:\n{e.output}", file=sys.stderr)
        return None
    line = None
    for ln in out.splitlines():
        if ln.startswith("MEASURE:"):
            line = ln
            break
    if line is None:
        print(f"[{backend}] runner output missing MEASURE: line.\nFull output:\n{out}", file=sys.stderr)
        return None
    bits = [int(tok) for tok in line.split(":", 1)[1].strip().split()]
    if len(bits) != meas_count:
        print(f"[{backend}] runner returned {len(bits)} bits but expected {meas_count}.", file=sys.stderr)
        return None
    return bits


def main():
    ap = argparse.ArgumentParser(description="Verify stab_cpu/stab_gpu measurement outcomes against Stim.")
    ap.add_argument("--n-qubits", type=int, required=True)
    ap.add_argument("--rounds", type=int, required=True)
    ap.add_argument("--ops-per-round", type=int, required=True)
    ap.add_argument("--seed", type=int, default=12345, help="Seed for circuit gen and Stim sampling.")
    ap.add_argument("--qasm-out", type=Path, default=Path("build/random_circuit.qasm"))
    ap.add_argument("--runner-exe", type=Path, default=Path("build/verify_from_qasm"))
    ap.add_argument("--backends", type=str, default="cpu,gpu", help="Comma list: cpu,gpu")
    args = ap.parse_args()

    args.qasm_out.parent.mkdir(parents=True, exist_ok=True)

    # 1) Build random circuit with Stim.
    stim_circ, stim_meas_count = build_random_stim_circuit(
        n_qubits=args.n_qubits,
        rounds=args.rounds,
        ops_per_round=args.ops_per_round,
        seed=args.seed,
    )

    # 2) Write OpenQASM 2.0 file.
    qasm_meas_count = write_openqasm(args.qasm_out, args.n_qubits, stim_circ)
    if qasm_meas_count != stim_meas_count:
        print("Internal error: measurement count mismatch between Stim and QASM export.", file=sys.stderr)
        sys.exit(2)

    # 3) Get Stim's one-shot measurement results (reproducible).
    sampler = stim_circ.compile_sampler(seed=args.seed)
    stim_bits = sampler.sample(1)[0].tolist()
    print(f"Stim produced {len(stim_bits)} measurement bits.")
    if len(stim_bits) != qasm_meas_count:
        print("Internal error: Stim sampler count mismatch.", file=sys.stderr)
        sys.exit(2)

    # 4) Optionally run CPU/GPU backends via provided runner.
    backends = [b.strip() for b in args.backends.split(",") if b.strip()]
    mismatches = 0
    for be in backends:
        if be not in ("cpu", "gpu", "nvgpu"):
            print(f"Skipping unknown backend '{be}'", file=sys.stderr)
            continue
        be_key = "nvgpu" if be == "gpu" else be
        cpp_bits = run_cpp_backend(args.runner_exe, be_key, args.qasm_out, args.seed, qasm_meas_count)
        if cpp_bits is None:
            continue
        # Compare bit-by-bit to Stim.
        diff_idx = [i for i, (a, b) in enumerate(zip(stim_bits, cpp_bits)) if a != b]
        if diff_idx:
            print(f"[{be_key}] {len(diff_idx)} mismatches vs Stim out of {qasm_meas_count}. Example first 10 idx: {diff_idx[:10]}")
            mismatches += len(diff_idx)
        else:
            print(f"[{be_key}] All measurements match Stim.")
    if mismatches > 0:
        print("Verification finished with mismatches. Note: RNG policies may differ between NWQ-Sim and Stim.")
        sys.exit(3)
    else:
        print("Verification passed.")
        sys.exit(0)


if __name__ == "__main__":
    main()
import argparse
import csv
import time
from pathlib import Path
import signal
import pyzx as zx

def count_t_gates(circ: zx.Circuit) -> int:
    # Prefer pyzx's built-in t-count if available
    try:
        return zx.tcount(circ)
    except Exception:
        cnt = 0
        for g in circ.gates:
            if getattr(g, "name", None) in ("T", "Tdg", "T*"):
                cnt += 1
        return cnt

def transpile_to_pbc(circ: zx.Circuit) -> zx.Circuit:
    # Work on a copy
    c = circ.copy()
    # Push to basic Clifford+T gate set
    try:
        c.to_basic_gates()
    except Exception:
        pass

    # Remove non-unitary gates that break to_graph on some PyZX versions
    try:
        c.gates = [
            g for g in c.gates
            if getattr(g, "name", "").upper() not in ("MEASURE", "BARRIER", "RESET")
        ]
    except Exception:
        pass

    # Convert via graph simplify + re-extract
    g = c.to_graph()
    zx.simplify.full_reduce(g)
    c = zx.extract.extract_circuit(g)
    return c

def transpile_to_pbc_with_timeout(circ: zx.Circuit, timeout_sec: int):
    # Linux/Unix-only timeout using SIGALRM
    def _handler(signum, frame):
        raise TimeoutError("Transpilation timed out")

    old_handler = signal.signal(signal.SIGALRM, _handler)
    try:
        signal.alarm(timeout_sec)
        t0 = time.perf_counter()
        pbc_circ = transpile_to_pbc(circ)
        dt_ms = (time.perf_counter() - t0) * 1000.0
        return pbc_circ, dt_ms
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

def process_folder(input_folder: Path, output_folder: Path, timeout_sec: int) -> None:
    output_folder.mkdir(parents=True, exist_ok=True)
    results_csv = output_folder / "results.csv"

    rows = [("circuit", "time_ms", "t_gates_after")]

    for qasm_path in sorted(input_folder.glob("*.qasm")):
        name = qasm_path.stem
        try:
            circ = zx.Circuit.from_qasm_file(str(qasm_path))
        except Exception as e:
            print(f"[WARN] Failed to parse {qasm_path.name}: {e}")
            continue

        try:
            pbc_circ, dt_ms = transpile_to_pbc_with_timeout(circ, timeout_sec)
        except TimeoutError:
            print(f"[WARN] {name}: transpilation exceeded {timeout_sec}s; skipping")
            rows.append((name, "TIMEOUT", ""))
            continue
        except KeyError as e:
            print(f"[WARN] {name}: to_graph failed ({e}); skipping")
            rows.append((name, "PARSE_ERROR", ""))
            continue

        t_after = count_t_gates(pbc_circ)

        # Save the transpiled (Pauli-based) circuit
        out_qasm = output_folder / f"{name}_pbc.qasm"
        try:
            with open(out_qasm, "w") as f:
                f.write(pbc_circ.to_qasm())
        except Exception as e:
            print(f"[WARN] Could not write {out_qasm.name}: {e}")

        rows.append((name, f"{dt_ms:.3f}", t_after))
        print(f"[OK] {name}: time={dt_ms:.2f} ms, T={t_after}")


def main():
    ap = argparse.ArgumentParser(description="Transpile Clifford+T QASM circuits to Pauli-based (PyZX t-par) and record stats.")
    ap.add_argument("--input", "-i", type=Path, required=True, help="Input folder containing .qasm files")
    ap.add_argument("--output", "-o", type=Path, required=True, help="Output folder for transpiled circuits and CSV")
    ap.add_argument("--timeout-minutes", type=int, default=10, help="Per-circuit transpilation timeout (minutes)")
    args = ap.parse_args()

    process_folder(args.input, args.output, args.timeout_minutes * 60)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt
import numpy as np


def load_stabsim_folder(folder: Path):
    """
    Parse files like multiplier_n45.txt where:
      - time_ms = (line1 + line2) * 1000
      - T_count = int(line4)
    Returns dict[name] = (time_ms, T_count)
    """
    out = {}
    for txt in sorted(folder.glob("*.txt")):
        try:
            lines = [x.strip() for x in txt.read_text().splitlines() if x.strip() != ""]
            if len(lines) < 4:
                print(f"[WARN] {txt.name}: not enough lines, skipping", file=sys.stderr)
                continue
            t_ms = (float(lines[0]) + float(lines[1])) * 1000.0
            t_count = int(lines[3])
            out[txt.stem] = (t_ms, t_count)
        except Exception as e:
            print(f"[WARN] Failed to parse {txt.name}: {e}", file=sys.stderr)
    return out


def load_results_csv(csv_path: Path):
    """
    Reads results.csv with columns: circuit,time_ms,t_gates_after
    Skips rows with TIMEOUT or missing values.
    Returns dict[name] = (time_ms, T_count)
    """
    out = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("circuit", "").strip()
            time_ms_s = (row.get("time_ms") or "").strip()
            t_s = (row.get("t_gates_after") or "").strip()
            if not name:
                continue
            if time_ms_s.upper() == "TIMEOUT":
                print(f"[INFO] {name}: timeout in results.csv, skipping", file=sys.stderr)
                continue
            try:
                time_ms = float(time_ms_s)
                t_count = int(t_s)
            except Exception:
                print(f"[WARN] {name}: bad row in results.csv (time_ms={time_ms_s}, t={t_s}), skipping", file=sys.stderr)
                continue
            out[name] = (time_ms, t_count)
    return out


def plot_comparison(pairs, out_path: Path, title: str):
    """
    pairs: list of (name, stabsim_ms, pyzx_ms)
    """
    names = [p[0] for p in pairs]
    stab_ms = [p[1] for p in pairs]
    pyzx_ms = [p[2] for p in pairs]

    n = len(names)
    if n == 0:
        print("[ERROR] No comparable circuits found (matching names and T counts).", file=sys.stderr)
        return

    fig_w = max(8.0, 0.4 * n + 2.0)
    plt.figure(figsize=(fig_w, 5.0))
    x = np.arange(n)
    width = 0.4

    plt.bar(x - width / 2, stab_ms, width, label="StabSim (ms)", color='red')
    plt.bar(x + width / 2, pyzx_ms, width, label="PyZX (ms)", color='blue')

    plt.xticks(x, names, rotation=45, ha="right")
    plt.ylabel("Time (ms)")
    plt.yscale('log')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"[OK] Saved plot to {out_path}")


def main():
    ap = argparse.ArgumentParser(description="Compare transpilation times between stabsim data and results.csv.")
    ap.add_argument("--stab-folder", type=Path, required=True,
                    help="Folder with *.txt files (e.g., stabsim_data/stab_T_time)")
    ap.add_argument("--results-csv", type=Path, required=True,
                    help="Path to results.csv produced by pyzx_transpiler")
    ap.add_argument("--out", type=Path, default=Path("time_comparison.png"),
                    help="Output image path (PNG)")
    args = ap.parse_args()

    stab = load_stabsim_folder(args.stab_folder)
    res = load_results_csv(args.results_csv)

    missing_in_res = sorted(set(stab.keys()) - set(res.keys()))
    if missing_in_res:
        print(f"[INFO] Missing in results.csv: {len(missing_in_res)} files (e.g., {missing_in_res[:5]})", file=sys.stderr)

    missing_in_stab = sorted(set(res.keys()) - set(stab.keys()))
    if missing_in_stab:
        print(f"[INFO] Missing in stabsim folder: {len(missing_in_stab)} files (e.g., {missing_in_stab[:5]})", file=sys.stderr)

    pairs = []
    mismatched_t = []
    for name in sorted(set(stab.keys()) & set(res.keys())):
        stab_ms, stab_t = stab[name]
        res_ms, res_t = res[name]
        if stab_t != res_t:
            mismatched_t.append((name, stab_t, res_t))
            continue
        pairs.append((name, stab_ms, res_ms))

    if mismatched_t:
        print(f"[WARN] T-count mismatch on {len(mismatched_t)} circuits; excluding from comparison.", file=sys.stderr)
        for name, s_t, r_t in mismatched_t[:10]:
            print(f"       {name}: stabsim T={s_t}, results.csv T={r_t}", file=sys.stderr)

    title = "Transpilation Time Comparison (ms) — matching T counts"
    plot_comparison(pairs, args.out, title)


if __name__ == "__main__":
    main()
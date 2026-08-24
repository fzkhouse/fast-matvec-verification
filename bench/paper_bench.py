#!/usr/bin/env python3
"""Run reference and original-baseline experiments on one paper parameter grid.

The data directory is explicit on purpose: large experiment fixtures are not
part of the release.  The command records only the verifier compute window;
file loading and table decoding occur before that window.
"""
import argparse
import csv
import re
import statistics
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def files_for(scheme: str, data: Path, rows: int):
    if scheme == "rpt":
        return (data / f"mat_B_{rows}.txt", None, data / f"vec_zu_{rows}.txt",
                data / f"mat_C_{rows}.txt", data / f"mat_ro_{rows}.txt", [])
    if scheme == "falcon512":
        return (data / f"Antt_{rows}.txt", data / f"Acoeff_{rows}.txt",
                data / f"z_{rows}.txt", data / f"u_{rows}.txt",
                data / f"table_{rows}.txt", ["--falcon-monty-table"])
    if scheme == "falcon1024":
        return (data / f"Antt_2_2_{rows}.txt", data / f"mat_A_coeff_2_2_{rows}.txt",
                data / f"z_2_2_{rows}.txt", data / f"u_2_2_{rows}.txt",
                data / f"table_2_2_1024_{rows}.txt", ["--falcon-monty-table"])
    if scheme == "kyber512":
        p = "5_9"
        return (data / f"mat_A_ntt_{p}_{rows}.txt", None,
                data / f"mat_z_{p}_{rows}.txt", data / f"mat_u_{p}_{rows}.txt",
                data / f"table_{p}_{rows}.txt", [])
    if scheme == "kyber1024":
        return (data / f"mat_A_ntt_4_{rows}.txt", None,
                data / f"mat_z_4_{rows}.txt", data / f"mat_u_4_{rows}.txt",
                data / f"table_4_{rows}.txt", [])
    if scheme in ("dilithium2", "dilithium5"):
        stem = f"{rows:04d}"
        return (data / f"mat_A_{stem}.txt", None,
                data / f"mat_z_{stem}.txt", data / f"mat_u_{stem}.txt",
                data / f"table_{stem}.txt", [])
    raise ValueError(scheme)


def command_for(scheme: str, implementation: str, files):
    a_ntt, a_coeff, z, u, table, extra = files
    if scheme == "rpt":
        return [ROOT / "build" / ("rpt-reference" if implementation == "reference" else "rpt-baseline"), a_ntt, z, u, table, "50"]
    if scheme.startswith("falcon"):
        if implementation == "reference":
            return [ROOT / "build/falcon-reference", scheme, a_ntt, z, u, table, *extra]
        return [ROOT / "build/falcon-baseline", a_coeff, z, u]
    if scheme.startswith("kyber"):
        binary = ROOT / "build" / ("kyber-reference" if implementation == "reference" else "kyber-baseline")
        if implementation == "reference":
            return [binary, scheme, a_ntt, z, u, table]
        return [binary, a_ntt, z, u]
    if scheme == "dilithium5":
        name = "dilithium-reference" if implementation == "reference" else "dilithium-baseline"
    else:
        name = f"{scheme}-reference" if implementation == "reference" else f"{scheme}-baseline"
    binary = ROOT / "build" / name
    if implementation == "reference":
        return [binary, z, u, table]
    return [binary, z, u, a_ntt, table]


def parse(stdout: str):
    passed = "result=PASS" in stdout or "Falcon-native matvec-eq: PASS" in stdout
    t = re.search(r"(?:compute_time|time_matvec_eq|time)=([0-9.]+)", stdout)
    c = re.search(r"(?:compute_cycles|cycles_matvec_eq|cycles)=([0-9]+)", stdout)
    kl = re.search(r"\bK=([0-9]+)\s+L=([0-9]+)", stdout)
    if not passed or not t:
        raise RuntimeError("benchmark did not pass or did not report metrics\n" + stdout)
    return float(t.group(1)), int(c.group(1)) if c else 0, int(kl.group(1)) if kl else 0, int(kl.group(2)) if kl else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme", required=True, choices=["falcon512", "falcon1024", "kyber512", "kyber1024", "dilithium2", "dilithium5", "rpt"])
    ap.add_argument("--rows", required=True, help="comma-separated row counts")
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--runs", type=int, default=9)
    ap.add_argument("--impl", choices=["reference", "baseline", "both"], default="both")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.runs < 1:
        ap.error("--runs must be positive")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows_list = [int(x) for x in args.rows.split(",")]
    implementations = ["reference", "baseline"] if args.impl == "both" else [args.impl]
    with args.out.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["scheme", "implementation", "rows", "K", "L", "runs", "pass_count", "median_time_s", "mean_time_s", "stddev_time_s", "median_cycles", "mean_cycles"])
        for rows in rows_list:
            files = files_for(args.scheme, args.data_root, rows)
            missing = [str(p) for p in files[:5] if p is not None and not p.is_file()]
            if missing:
                raise FileNotFoundError("missing fixture(s): " + ", ".join(missing))
            for implementation in implementations:
                times, cycles = [], []
                K = L = 0
                for _ in range(args.runs):
                    cp = subprocess.run([str(x) for x in command_for(args.scheme, implementation, files)], text=True, capture_output=True)
                    if cp.returncode != 0:
                        raise RuntimeError(cp.stdout + cp.stderr)
                    t, c, K, L = parse(cp.stdout)
                    times.append(t)
                    cycles.append(c)
                writer.writerow([args.scheme, implementation, rows, K, L, args.runs, len(times),
                                 f"{statistics.median(times):.6f}", f"{statistics.mean(times):.6f}",
                                 f"{statistics.pstdev(times):.6f}", int(statistics.median(cycles)),
                                 f"{statistics.mean(cycles):.0f}"])
                out.flush()


if __name__ == "__main__":
    main()

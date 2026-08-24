#!/usr/bin/env python3
"""Generate deterministic matrix-vector relation fixtures.

The text format is deliberately simple: one header ``q n rows cols`` followed
by one polynomial per line.  A is K x L, z is rows x L, and u is rows x K.
"""
import argparse
import os
import random

PARAMS = {
    "falcon512": (12289, 512, 1, 1),
    "falcon1024": (12289, 1024, 2, 2),
    "kyber512": (3329, 256, 2, 2),
    "kyber1024": (3329, 256, 4, 4),
    "dilithium2": (8380417, 256, 4, 4),
    "dilithium5": (8380417, 256, 8, 7),
}

def add(a, b, q):
    return [(x + y) % q for x, y in zip(a, b)]

def mul(a, b, q):
    n = len(a)
    out = [0] * n
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                p = i + j
                out[p % n] = (out[p % n] + x * y * (1 if p < n else -1)) % q
    return out

def write_grid(path, q, n, rows, cols, polys):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{q} {n} {rows} {cols}\n")
        for p in polys:
            f.write(" ".join(map(str, p)) + "\n")

def generate(args):
    q, n, k, l = PARAMS[args.scheme]
    if args.rows < 1:
        raise SystemExit("rows must be positive")
    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)

    # The default sparse public fixture keeps large experiments practical
    # without external numerical packages. --dense opts into schoolbook
    # generation for small, more demanding fixtures.
    A = []
    for i in range(k):
        for j in range(l):
            p = [0] * n
            if not args.dense and i == j:
                p[0] = 1
            else:
                p = [rng.randrange(q) for _ in range(n)]
            A.append(p)

    z = [[rng.randrange(q) for _ in range(n)] for _ in range(args.rows * l)]
    u = []
    for r in range(args.rows):
        for i in range(k):
            acc = [0] * n
            for j in range(l):
                acc = add(acc, mul(A[i * l + j], z[r * l + j], q), q)
            u.append(acc)

    write_grid(os.path.join(args.out, "A.txt"), q, n, k, l, A)
    write_grid(os.path.join(args.out, "z.txt"), q, n, args.rows, l, z)
    write_grid(os.path.join(args.out, "u.txt"), q, n, args.rows, k, u)
    with open(os.path.join(args.out, "meta.txt"), "w", encoding="utf-8") as f:
        f.write(f"scheme={args.scheme}\nrows={args.rows}\nseed={args.seed}\n")
    print(f"generated scheme={args.scheme} rows={args.rows} seed={args.seed} out={args.out}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scheme", choices=sorted(PARAMS), required=True)
    p.add_argument("--rows", type=int, required=True)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--out", required=True)
    p.add_argument("--dense", action="store_true")
    generate(p.parse_args())

if __name__ == "__main__":
    main()

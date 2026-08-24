#!/usr/bin/env python3
import os
import math
import argparse
import random
from typing import List, Tuple, Optional

# ---------------- basic math ----------------

def mod_pow(a: int, e: int, q: int) -> int:
    return pow(a, e, q)

def factorize(x: int) -> List[int]:
    f = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            f.append(p)
            while x % p == 0:
                x //= p
        p += 1
    if x > 1:
        f.append(x)
    return f

def primitive_root(q: int) -> int:
    phi = q - 1
    fac = factorize(phi)
    for g in range(2, q):
        ok = True
        for p in fac:
            if mod_pow(g, phi // p, q) == 1:
                ok = False
                break
        if ok:
            return g
    raise ValueError("no primitive root found")

def bit_reverse(a: List[int]) -> None:
    n = len(a)
    j = 0
    i = 1
    while i < n:
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
        i += 1

def ntt_ctx_init(n: int, q: int):
    if n <= 0 or (n & (n - 1)) != 0:
        raise ValueError("n must be power of two")
    if (q - 1) % (2 * n) != 0:
        raise ValueError(f"(q-1)%(2n)!=0, this NTT impl not support q={q}, n={n}")

    g = primitive_root(q)
    exp = (q - 1) // (2 * n)
    psi = mod_pow(g, exp, q)
    omega = (psi * psi) % q
    if mod_pow(psi, n, q) != (q - 1) % q:
        raise ValueError("psi^n != -1 mod q")
    if mod_pow(omega, n, q) != 1:
        raise ValueError("omega^n != 1 mod q")

    psi_pows = [1] * n
    for i in range(1, n):
        psi_pows[i] = (psi_pows[i - 1] * psi) % q
    return {"n": n, "q": q, "omega": omega, "psi_pows": psi_pows}

def ntt_forward(poly: List[int], ctx) -> List[int]:
    # ULA-NTT layout: DIF forward, bit-reversed output (no explicit bit_reverse)
    n = ctx["n"]
    q = ctx["q"]
    omega = ctx["omega"]
    psi_pows = ctx["psi_pows"]

    a = [0] * n
    for i in range(n):
        a[i] = (poly[i] % q) * psi_pows[i] % q

    length = n
    while length >= 2:
        step = n // length
        wlen = mod_pow(omega, step, q)
        half = length >> 1
        for i in range(0, n, length):
            w = 1
            for j in range(half):
                u = a[i + j]
                v = a[i + j + half]
                x = u + v
                if x >= q:
                    x -= q
                y = u - v
                if y < 0:
                    y += q
                a[i + j] = x
                a[i + j + half] = (y * w) % q
                w = (w * wlen) % q
        length >>= 1
    return a

# ---------------- Falcon Montgomery (q=12289) ----------------

FALCON_Q = 12289
FALCON_Q0I = 12287   # -1/q mod 2^16
FALCON_R2 = 10952    # 2^32 mod q

def falcon_monty_mul_12289(x: int, y: int) -> int:
    z = (x * y) & 0xFFFFFFFF
    w = (((z * FALCON_Q0I) & 0xFFFF) * FALCON_Q) & 0xFFFFFFFF
    z = (z + w) >> 16
    z -= FALCON_Q
    if z < 0:
        z += FALCON_Q
    return z

def falcon_to_monty_12289(x: int) -> int:
    return falcon_monty_mul_12289(x % FALCON_Q, FALCON_R2)

def vec_to_monty_12289(v: List[int]) -> List[int]:
    return [falcon_to_monty_12289(x) for x in v]

# ---------------- IO format ----------------
# grid format:
# q n rows cols
# then rows*cols lines, each line is one polynomial of length n

def write_grid(path: str, q: int, n: int, rows: int, cols: int, polys: List[List[int]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{q} {n} {rows} {cols}\n")
        idx = 0
        for _ in range(rows):
            for _ in range(cols):
                f.write(" ".join(map(str, polys[idx])) + "\n")
                idx += 1

def read_grid(path: str) -> Tuple[int, int, int, int, List[List[int]]]:
    with open(path, "r", encoding="utf-8") as f:
        head = f.readline().strip().split()
        if len(head) != 4:
            raise ValueError(f"bad grid header: {path}")
        q, n, rows, cols = map(int, head)
        polys = []
        for _ in range(rows * cols):
            vals = f.readline().strip().split()
            if len(vals) != n:
                raise ValueError(f"bad poly line length in {path}")
            polys.append([int(x) % q for x in vals])
    return q, n, rows, cols, polys

# table format:
# q n R T K L
# for each row r:
#   randv[0..T-1]  (one line)
#   gamma_ntt[k] for k in [0,K)  (K lines, each length n)
#   AT_gamma_ntt[l] for l in [0,L) (L lines, each length n)

def read_table_header_body(path: str):
    if not os.path.exists(path):
        return None, ""
    with open(path, "r", encoding="utf-8") as f:
        first = f.readline().strip().split()
        if len(first) != 6:
            return "BAD_HEADER", ""
        q, n, R, T, K, L = map(int, first)
        body = f.read()
    return (q, n, R, T, K, L), body

# ---------------- dataset generation ----------------

def params_from_scheme(name: str):
    if name == "falcon512":
        return {"q": 12289, "n": 512, "K": 1, "L": 1}
    if name == "falcon1024":
        return {"q": 12289, "n": 512, "K": 2, "L": 2}
    # 当前 NTT 实现要求 (q-1)%(2n)==0，下面两组会被拦截
    if name == "kyber":
        return {"q": 3329, "n": 256, "K": 2, "L": 2}
    if name == "dilithium":
        return {"q": 8380417, "n": 256, "K": 4, "L": 4}
    raise ValueError(f"unknown scheme: {name}")

def rand_poly(rng: random.Random, n: int, q: int) -> List[int]:
    return [rng.randrange(q) for _ in range(n)]

def mul_negacyclic(a: List[int], b: List[int], q: int) -> List[int]:
    n = len(a)
    out = [0] * n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            t = ai * bj
            s = i + j
            if s < n:
                out[s] += t
            else:
                out[s - n] -= t
    return [x % q for x in out]

def add_poly(a: List[int], b: List[int], q: int) -> List[int]:
    return [(x + y) % q for x, y in zip(a, b)]

def is_falcon_scheme(scheme: str) -> bool:
    return scheme in ("falcon512", "falcon1024")

def resolve_bool_flag(enable_opt: bool, disable_opt: bool, default_val: bool) -> bool:
    if enable_opt and disable_opt:
        raise ValueError("conflicting flags: both enable and disable are set")
    if enable_opt:
        return True
    if disable_opt:
        return False
    return default_val

def _resolve_out_path(out_dir: str, name_or_path: str) -> str:
    # 绝对路径就直接用；相对路径则拼到 --out 目录下
    if os.path.isabs(name_or_path):
        return name_or_path
    return os.path.join(out_dir, name_or_path)

def gen_azu(args):
    P = params_from_scheme(args.scheme)
    q, n = P["q"], P["n"]
    K = args.k if args.k is not None else P["K"]
    L = args.l if args.l is not None else P["L"]
    N_rows = args.nrows

    if not (1 <= K <= 10 and 1 <= L <= 10):
        raise ValueError("K,L must be in [1,10]")

    use_falcon_monty_a = resolve_bool_flag(
        args.falcon_monty_a_ntt, args.no_falcon_monty_a_ntt, is_falcon_scheme(args.scheme)
    )
    if use_falcon_monty_a and q != FALCON_Q:
        raise ValueError("Falcon monty A_ntt mode requires q=12289")

    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)

    ctx = ntt_ctx_init(n, q)

    # A in coeff for building u
    A_coeff = [rand_poly(rng, n, q) for _ in range(K * L)]

    # save A in NTT or NTT+Monty
    A_ntt_raw = [ntt_forward(p, ctx) for p in A_coeff]
    if use_falcon_monty_a:
        A_ntt = [vec_to_monty_12289(p) for p in A_ntt_raw]
    else:
        A_ntt = A_ntt_raw

    # z in coeff
    z = [rand_poly(rng, n, q) for _ in range(N_rows * L)]

    # u = A * z in coeff domain (always from A_coeff)
    u = []
    for i in range(N_rows):
        for k in range(K):
            acc = [0] * n
            for l in range(L):
                prod = mul_negacyclic(A_coeff[k * L + l], z[i * L + l], q)
                acc = add_poly(acc, prod, q)
            u.append(acc)

    # 可自定义输出文件名
    path_a_coeff = _resolve_out_path(args.out, args.file_a_coeff)
    path_a_ntt   = _resolve_out_path(args.out, args.file_a_ntt)
    path_z       = _resolve_out_path(args.out, args.file_z)
    path_u       = _resolve_out_path(args.out, args.file_u)

    write_grid(path_a_coeff, q, n, K, L, A_coeff)
    write_grid(path_a_ntt, q, n, K, L, A_ntt)
    write_grid(path_z, q, n, N_rows, L, z)
    write_grid(path_u, q, n, N_rows, K, u)

    print("generated:")
    print(path_a_coeff)
    print(path_a_ntt)
    print(path_z)
    print(path_u)
    print(f"[mode] A_ntt_falcon_monty={'ON' if use_falcon_monty_a else 'OFF'}")

def gen_table(args):
    P = params_from_scheme(args.scheme)
    q_s, n_s = P["q"], P["n"]

    q, n, K, L, A_ntt = read_grid(args.mat_a_ntt)
    if q != q_s or n != n_s:
        raise ValueError("scheme and mat_A_ntt header mismatch")

    ctx = ntt_ctx_init(n, q)  # for gamma coeff->NTT
    rng = random.Random(args.seed)

    use_falcon_monty_table = resolve_bool_flag(
        args.falcon_monty_table, args.no_falcon_monty_table, is_falcon_scheme(args.scheme)
    )
    if use_falcon_monty_table and q != FALCON_Q:
        raise ValueError("Falcon monty table mode requires q=12289")

    A_is_monty = resolve_bool_flag(
        args.a_ntt_is_monty, args.a_ntt_is_non_monty, use_falcon_monty_table
    )

    if args.rows_add is None or args.rows_add <= 0:
        raise ValueError("--rows-add must be explicitly set to a positive integer")
    rows_add = args.rows_add

    old = read_table_header_body(args.table)
    if old[0] == "BAD_HEADER":
        if not args.reset_table:
            raise ValueError(f"bad table header in {args.table}; use --reset-table")
        old_body = ""
        old_R = 0
        T = args.table_t
    elif old[0] is None:
        old_body = ""
        old_R = 0
        T = args.table_t
    else:
        (oq, on, oR, oT, oK, oL), old_body = old
        T = oT
        if (oq, on, oK, oL) != (q, n, K, L):
            if not args.reset_table:
                raise ValueError("existing table header mismatch; use --reset-table")
            old_body = ""
            old_R = 0
            T = args.table_t
        else:
            old_R = oR

    if T <= 0:
        raise ValueError("table T must be positive")

    if old_R > 0 and not args.reset_table:
        print("[warn] appending to existing table; ensure representation mode is same as old content.")

    new_R = old_R + rows_add
    lines = []

    for _ in range(rows_add):
        randv = [rng.randrange(q) for _ in range(T)]
        lines.append(" ".join(map(str, randv)))

        gamma_ntt_raw = []
        gamma_ntt_out = []

        for _k in range(K):
            g_coeff = rand_poly(rng, n, q)
            g_ntt = ntt_forward(g_coeff, ctx)     # normal NTT
            gamma_ntt_raw.append(g_ntt)

            if use_falcon_monty_table:
                g_store = vec_to_monty_12289(g_ntt)   # NTT+Monty
            else:
                g_store = g_ntt
            gamma_ntt_out.append(g_store)
            lines.append(" ".join(map(str, g_store)))

        # Build AT_gamma_ntt[l]
        for l in range(L):
            if use_falcon_monty_table:
                if A_is_monty:
                    # A_monty * gamma_monty --montyMul--> monty ; sum keeps monty
                    acc_monty = [0] * n
                    for k in range(K):
                        Akl = A_ntt[k * L + l]           # monty
                        gk_m = gamma_ntt_out[k]          # monty
                        for i in range(n):
                            prod_m = falcon_monty_mul_12289(Akl[i], gk_m[i])  # monty
                            acc_monty[i] += prod_m
                            if acc_monty[i] >= q:
                                acc_monty[i] -= q
                    lines.append(" ".join(map(str, acc_monty)))
                else:
                    # A_normal * gamma_normal -> normal ; then to monty
                    acc_norm = [0] * n
                    for k in range(K):
                        Akl = A_ntt[k * L + l]           # normal
                        gk_n = gamma_ntt_raw[k]          # normal
                        for i in range(n):
                            acc_norm[i] = (acc_norm[i] + Akl[i] * gk_n[i]) % q
                    acc_monty = vec_to_monty_12289(acc_norm)
                    lines.append(" ".join(map(str, acc_monty)))
            else:
                # normal mode
                acc = [0] * n
                for k in range(K):
                    Akl = A_ntt[k * L + l]
                    gk = gamma_ntt_out[k]
                    for i in range(n):
                        acc[i] = (acc[i] + Akl[i] * gk[i]) % q
                lines.append(" ".join(map(str, acc)))

    os.makedirs(os.path.dirname(args.table) or ".", exist_ok=True)
    with open(args.table, "w", encoding="utf-8") as f:
        f.write(f"{q} {n} {new_R} {T} {K} {L}\n")
        if old_body:
            f.write(old_body)
            if not old_body.endswith("\n"):
                f.write("\n")
        f.write("\n".join(lines))
        if lines:
            f.write("\n")

    print("table updated:", args.table)
    print(f"rows old={old_R}, add={rows_add}, now={new_R}")
    print(f"[mode] falcon_monty_table={'ON' if use_falcon_monty_table else 'OFF'}")
    if use_falcon_monty_table:
        print(f"[mode] input_A_ntt_is_monty={'ON' if A_is_monty else 'OFF'}")

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("gen-azu", help="generate A_ntt, z, u")
    p1.add_argument("--scheme", default="falcon512", choices=["falcon512", "falcon1024", "kyber", "dilithium"])
    p1.add_argument("--nrows", type=int, default=64)
    p1.add_argument("--k", type=int, default=None)
    p1.add_argument("--l", type=int, default=None)
    p1.add_argument("--seed", type=int, default=12345)
    p1.add_argument("--out", default="data")
    p1.add_argument("--falcon-monty-a-ntt", action="store_true", help="save A_ntt in Falcon NTT+Monty")
    p1.add_argument("--no-falcon-monty-a-ntt", action="store_true", help="force save A_ntt in normal NTT")

    # 新增：自定义输出文件名（默认兼容旧逻辑）
    p1.add_argument("--file-a-coeff", default="mat_A_coeff.txt")
    p1.add_argument("--file-a-ntt",   default="mat_A_ntt.txt")
    p1.add_argument("--file-z",       default="mat_z.txt")
    p1.add_argument("--file-u",       default="mat_u.txt")

    p1.set_defaults(func=gen_azu)

    p2 = sub.add_parser("gen-table", help="append-only table generation from A_ntt")
    p2.add_argument("--scheme", default="falcon512", choices=["falcon512", "falcon1024", "kyber", "dilithium"])
    p2.add_argument("--mat-a-ntt", required=True)
    p2.add_argument("--table", required=True)
    p2.add_argument("--rows-add", type=int, required=True, help="number of rows to append")
    p2.add_argument("--table-t", type=int, default=1, help="randv length T for NEW table (ignored when appending)")
    p2.add_argument("--seed", type=int, default=2026)
    p2.add_argument("--reset-table", action="store_true", help="ignore existing table content and start fresh")
    p2.add_argument("--falcon-monty-table", action="store_true",
                    help="generate gamma_ntt and AT_gamma_ntt in Falcon NTT+Monty")
    p2.add_argument("--no-falcon-monty-table", action="store_true",
                    help="force normal (non-Monty) table")

    p2.add_argument("--a-ntt-is-monty", action="store_true",
                    help="tell generator input mat_A_ntt is already Falcon NTT+Monty")
    p2.add_argument("--a-ntt-is-non-monty", action="store_true",
                    help="tell generator input mat_A_ntt is normal NTT")
    p2.set_defaults(func=gen_table)

    args = ap.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
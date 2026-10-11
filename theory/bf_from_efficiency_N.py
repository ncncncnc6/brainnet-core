#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/bf_from_efficiency_N.py —— P1′：N 靶点树模型下由 ε 反推 bf（theory-p1b/1.0.3）

P1 的否证：2 靶点玩具的 ε 下界 = 0.9017，够不到实测 ε_spiny=0.86、ε_basket=0.76。
P1′ 的假设：**靶点数越多、可共享的主干越深，线长效率的下界越低**，于是实测值可达。
本件把 N 从 2 扫到 64 检验该假设，并在可达的 N 上反推 bf。

模型（单参数 t，定义清楚）
  根在原点；N 个靶点均匀分布在半径 R、半角 ±Φ 的扇面上。
  主干沿角平分线到 t·R；其后每个靶点各接一支：b_i(t) = |p_i − trunk_end|
    wire(t) = t·R + Σ_i b_i(t)
    path(t) = Σ_i (t·R + b_i(t)) = N·t·R + Σ_i b_i(t)      ← 主干被 N 条路径各走一次
  极小化 C(t) = wire(t) + bf · path(t)；效率 ε = wire_min/wire(t*)、γ = path_min/path(t*)

判据
  [1] ε 下界随 N 单调下降（假设成立？）
  [2] 存在最小 N*，使 ε=0.86（spiny）与 0.76（basket）可达；在 N* 上给出对应 bf
  [3] 平衡点（ε=γ）处的 bf 随 N 的变化
  [4] 负例：N=1（无共享）→ 无权衡（ε 与 γ 都不随 bf 变）
"""
VERSION = "theory-p1b/1.0.3"
import math, sys

R, PHI = 1.0, 0.6

def targets(N, phi=PHI):
    if N == 1: return [(R, 0.0)]
    return [(R * math.cos(-phi + 2 * phi * i / (N - 1)),
             R * math.sin(-phi + 2 * phi * i / (N - 1))) for i in range(N)]

def wire_path(t, pts):
    bx = t * R                                     # 主干末端（角平分线 = x 轴）
    trunk = t * R
    brs = [math.hypot(p[0] - bx, p[1]) for p in pts]
    wire = trunk + sum(brs)
    path = len(pts) * trunk + sum(brs)
    return wire, path

def sweep(N, samples=800):
    pts = targets(N)
    ws = [wire_path(i / samples, pts) for i in range(samples + 1)]
    wmin = min(w for w, _ in ws); pmin = min(p for _, p in ws)
    return pts, ws, wmin, pmin

def efficiencies(bf, pts, ws, wmin, pmin, samples=800):
    best = min(range(len(ws)), key=lambda i: ws[i][0] + bf * ws[i][1])
    w, p = ws[best]
    return wmin / w, pmin / p, best / samples

def eps_floor(N):
    """bf→∞ 时的 ε（完全按路径优先 → t=0 或 t=1 端），以及 bf=0 时的 ε。"""
    pts, ws, wmin, pmin = sweep(N)
    e_inf, _, _ = efficiencies(1e6, pts, ws, wmin, pmin)
    e_0, _, _ = efficiencies(0.0, pts, ws, wmin, pmin)
    return min(e_0, e_inf), e_0, e_inf

def bf_for_eps(target, N, lo=0.0, hi=1e6, iters=70):
    pts, ws, wmin, pmin = sweep(N)
    for _ in range(iters):
        mid = (lo + hi) / 2
        e, _, _ = efficiencies(mid, pts, ws, wmin, pmin)
        if e > target: lo = mid
        else: hi = mid
    bf = (lo + hi) / 2
    return bf, efficiencies(bf, pts, ws, wmin, pmin)[0]

def bf_for_eps_bounded(target, N, hi=1.0, iters=60):
    """[D 组 (b) 5-6/58-66/92-102] **把文献定义域 bf ∈ [0,1] 当硬约束**：
    只在 [0,1] 内二分；若连 bf=1 都够不到目标，返回 None（不可达），不再给出定义域外的 bf。
    理由：PMC11790248 逐字定义 bf=0 只最小化电缆、bf=1 电缆与延迟等权；越域值不是同一个参数。
    """
    pts, ws, wmin, pmin = sweep(N)
    e1, _, _ = efficiencies(hi, pts, ws, wmin, pmin)
    if e1 > target: return None                      # 域内不可达
    lo = 0.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        e, _, _ = efficiencies(mid, pts, ws, wmin, pmin)
        if e > target: lo = mid
        else: hi = mid
    bf = (lo + hi) / 2
    return bf, efficiencies(bf, pts, ws, wmin, pmin)[0]

def balance_bf(N, lo=0.0, hi=1e6, iters=70):
    pts, ws, wmin, pmin = sweep(N)
    for _ in range(iters):
        mid = (lo + hi) / 2
        e, g, _ = efficiencies(mid, pts, ws, wmin, pmin)
        if e > g: lo = mid
        else: hi = mid
    bf = (lo + hi) / 2
    e, g, _ = efficiencies(bf, pts, ws, wmin, pmin)
    return bf, e, g

def selftest():
    print(f"{VERSION}\n")
    ok = True
    print("[1] ε 下界随 N 的变化（假设：N 越大下界越低）")
    floors = []
    for N in (2, 4, 8, 16, 32, 64):
        fl, e0, einf = eps_floor(N)
        floors.append((N, fl))
        print(f"    N={N:<3} ε(bf=0)={e0:.4f}  ε(bf→∞)={einf:.4f}  → 下界={fl:.4f}")
    mono = all(floors[i][1] >= floors[i + 1][1] - 1e-9 for i in range(len(floors) - 1))
    print(f"    下界随 N 单调下降: {mono}")
    ok &= mono

    print("\n[2] 实测 ε 的**域内**可达性与反推 bf（bf 限制在文献定义域 [0,1]）")
    spiny_ok = False
    for tag, target in (("ε_spiny", 0.86), ("ε_basket", 0.76), ("ε_star", 0.02)):
        hit = None
        for N in (2, 4, 8, 16, 32, 64, 128):
            r = bf_for_eps_bounded(target, N)
            if r is not None: hit = (N, r[0], r[1]); break
        if hit is None:
            fl = eps_floor(128)[0]
            print(f"    {tag}={target}：在 bf≤1 的约束下扫到 N=128 仍不可达（N=128 的 ε 下界 {fl:.4f}）→ 如实报告")
        else:
            N_, bf_, e_ = hit
            print(f"    {tag}={target}：域内最小可达 N*={N_} ⇒ bf ≈ {bf_:.3f}（≤1，复核 ε={e_:.4f}）")
            if tag == "ε_spiny": spiny_ok = True
    print(f"    域内结论：ε_spiny 在 bf≤1 下可达 = {spiny_ok}（ε_basket 需越域，见上）")
    ok &= spiny_ok      # [D 组 (b) 92-108] 原版 [2] 只打印、不进总判 → 现并入

    print("\n[3] 平衡点（ε=γ）处的 bf 随 N（**并入总判**）")
    in_domain = True
    for N in (2, 8, 32):
        bf, e, g = balance_bf(N)
        d_ok = (0.0 <= bf <= 1.0)
        in_domain &= d_ok
        print(f"    N={N:<3} bf≈{bf:.3f}（ε={e:.4f}, γ={g:.4f}）  落在定义域 [0,1]: {d_ok}")
    print(f"    平衡点全部落在文献定义域 [0,1] 内: {in_domain}")
    ok &= in_domain

    print("\n[4] 负例：N=1（无共享）→ 无权衡")
    pts, ws, wmin, pmin = sweep(1)
    es = [efficiencies(bf, pts, ws, wmin, pmin)[0] for bf in (0.0, 1.0, 10.0, 1e3)]
    flat = max(es) - min(es) < 1e-9
    caught = flat; ok &= caught
    print(f"    N=1 时 ε = {[round(x,4) for x in es]} → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

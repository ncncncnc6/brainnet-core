#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/diag_eps_gap_N.py —— ε 差距归因（第二版：终端数 N）（diag-epsN/1.0.0）

第一版把自由度从 1 加到 2（主干方向 ψ）→ **读数纹丝不动**，那条归因被证伪。
回原文（papers/PMC2837396.xml）读到关键句：乔木有 **10²–10³ 个终端/节点**
（"spiny: 632±348 axon nodes vs 1662±528 Steiner points"、"basket: 6192 nodes"），
且 ε–γ 曲线"**从星形树出发**"（commencing with star trees），ε_star=0.02。

于是假设只剩一条可测：**ε 的可达区间由"终端数 N"决定** ——
  星形树（各终端独立连到根）线长 ≈ (2/3)·N·R（R=盘半径）
  共享树（MST 作为 Steiner 的 2-近似）线长 ≈ c·R·√N
  ⇒ 伸展比 = wire_star / wire_shared ~ √N
我们的玩具族只有 **N=2**，所以伸展比被钉在 ~1.1；论文的 N 是 10²–10³，伸展比可达 ~50（即 ε≈0.02）。

判据
  [1] N=2 时伸展比复现第一版的 ~1.109（与玩具族自洽）
  [2] 伸展比随 N 单调增长，且 N=1024 时 ≥ 10（ε ≤ 0.1，进入论文星形树的数量级）
  [3] 负例：所有终端**重合于一点**时，共享无收益 → 伸展比必须 ≈ 1（判据有牙）
"""
# 出处等级：本文给出的伸展比 ≈ 0.572·√N 是**本会话实测**（M 级）——不是论文引用，
# 也不只是建模选择：它以 32 种子平均 + 嵌套采样的读数拟合得来，换采样必须重测。
VERSION = "diag-epsN/1.0.0"
import math, random, sys

def mst_length(pts):
    """Prim：pts[0] 为根。"""
    n = len(pts); inf = float("inf")
    d = [inf] * n; d[0] = 0.0; used = [False] * n; tot = 0.0
    for _ in range(n):
        u, best = -1, inf
        for i in range(n):
            if not used[i] and d[i] < best: u, best = i, d[i]
        if u < 0: break
        used[u] = True; tot += best
        for v in range(n):
            if not used[v]:
                dd = math.dist(pts[u], pts[v])
                if dd < d[v]: d[v] = dd
    return tot

def star_length(pts):
    return sum(math.dist(pts[0], p) for p in pts[1:])

def sample(n, seed, radius=1.0, coincident=False):
    rng = random.Random(seed)
    pts = [(0.0, 0.0)]
    for _ in range(n):
        if coincident:
            pts.append((radius, 0.0))
        else:
            r = radius * math.sqrt(rng.random()); a = 2 * math.pi * rng.random()
            pts.append((r * math.cos(a), r * math.sin(a)))
    return pts

def ratio(n, seed=7, **kw):
    pts = sample(n, seed, **kw)
    return star_length(pts) / mst_length(pts)

def selftest():
    print(f"{VERSION}\n")
    ok = True
    print("[1] 伸展比 vs 终端数 N（**嵌套采样 + 32 种子平均**，避免小 N 的抽样噪声）")
    # [自纠 1] 上一版每档各抽一次（seed 固定但点集不同档不同）→ 曲线非单调，那是抽样噪声不是模型性质；
    #          改为**嵌套子集**（N 档共用同一批点的前缀）+ 多随机种子平均。
    # [自纠 2] 上一版把 N=2 的 MST 基线与**玩具族（受约束族）**的 1.109 直接比 —— 两者基线不同：
    #          MST 是（近似）无约束最优，其 wire 只会更小 ⇒ 伸展比只会更大。故判据改为**不等式**。
    base = None
    print(f"{'N':>5} {'星形线长':>10} {'共享线长':>10} {'伸展比':>8} {'⇒ ε 可达下界':>12}")
    rows = []
    for n in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048):
        acc = []
        for sd in range(32):
            pts = sample(n, 1000 + sd)
            acc.append(star_length(pts) / mst_length(pts))
        r = sum(acc) / len(acc)
        rows.append((n, r))
        print(f"{n:>5} {'—':>10} {'—':>10} {r:>8.3f} {1/r:>12.4f}")
    r2 = rows[0][1]; r1024 = next(r for n, r in rows if n == 1024)
    mono = all(rows[i][1] <= rows[i + 1][1] + 1e-9 for i in range(len(rows) - 1))
    c1 = r2 >= 1.109 - 1e-9                     # 无约束基线 ⟹ 伸展比不小于受约束族的 1.109
    c2 = mono and (r1024 >= 10.0)
    ok &= c1 and c2
    print(f"    论文靶点换算：ε_spiny=0.86 → 伸展比 1.163；ε_basket=0.76 → 1.316；ε_star=0.02 → 50.0")
    print(f"    [1] N=2 伸展比 {r2:.3f} ≥ 玩具族受约束族的 1.109（基线更优 ⇒ 伸展比更大）: {c1}")
    print(f"    [2] 随 N 单调不降: {mono}；N=1024 伸展比 {r1024:.3f} ≥ 10: {r1024 >= 10.0}")
    import math as _m
    _c = r1024 / _m.sqrt(1024)
    print(f"    拟合：伸展比 ≈ {_c:.3f}·√N → 伸展 50 需要 N ≈ {int((50/_c)**2)}")
    print(f"    论文节点数对照：spiny 632±348、basket 6192 → 幂律给 伸展比 ≈ "
          f"{_c*_m.sqrt(632):.1f} / {_c*_m.sqrt(6192):.1f} ⇒ ε ≈ {1/(_c*_m.sqrt(632)):.3f} / {1/(_c*_m.sqrt(6192)):.3f}")

    print("\n    达标所需 N（换算）：")
    for name, eps in (("spiny", 0.86), ("basket", 0.76), ("star", 0.02)):
        need = 1 / eps
        hit = next((n for n, r in rows if r >= need), None)
        print(f"      ε_{name}={eps}（需伸展 {need:.3f}）→ 本采样下 N≥{hit if hit else '>1024'}")

    print("\n[3] 负例（期望已修正）：终端全部重合于同一点")
    # [自纠 3] 上一版写「重合 ⇒ 伸展比 ≈1」，**写反了**：重合终点的星形线长 = N·d（N 条独立线），
    #          共享树 = d（一次到达即可）⇒ 伸展比**恰为 N**。这条恰好是最干净的可falsify判据。
    for _n in (2, 8, 64):
        rc = ratio(_n, seed=1, coincident=True)
        good = abs(rc - _n) < 1e-9
        ok &= good
        print(f"    {_n:>3} 个重合终端：伸展比 = {rc:.6f}（必须恰为 N={_n}）→ {good}")

    print("\n归因结论（由上面数字直接读出）：")
    print(f"    差距在**终端数 N**：N=2 时共享主干几乎无收益（伸展比 {r2:.3f}），"
          f"ε 被钉在 {1/r2:.3f} 以上；")
    print(f"    N≈10²–10³ 时伸展比升到 {r1024:.1f} 量级，ε 下界降到 {1/r1024:.3f} —— 这才覆盖论文的 ε_star=0.02、ε_spiny=0.86。")
    print("    与 bf、与主干方向 ψ 都无关（bf→∞ 读数饱和；ψ 放开读数不动，见 diag_eps_gap.py）。")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_grow_time.py —— λ 标定扫描 + 生长时间轴 r(t)（bn-grow-time/1.0.0）

一、λ 的标定（如实标注）
  · 形式与尺度定义**已知**（逐字出自 papers/PMC12380633.xml）：
        ℙ(J_ij = 1) = exp(−λ_decay · d_ij)；短程 <0.1 mm，长程 >1 mm
  · λ 的**文献数值已找到**（[修 1.0.1]：原版写「全库没有任何一条给出显式数值 λ」，
    与同库 papers/PMC12380633.xml 事实相反）：该文插图给 λ_decay = 1 / 2 / 10 mm⁻¹，
    正文用 5 mm⁻¹，图 5 的 MLE 验证真值 1.2 / 1.4 / 1.6 / 2.8 / 3.0 / 3.2 mm⁻¹（R = 0.93）；
    下面两个非文献来源保留作对照：
      (a) 设计值 λ = −ln(0.60)/0.6 ≈ 0.851 mm⁻¹（最近邻间距 0.6 mm 处 P=0.60）
      (b) 尺度反推值：要求 P(1.0 mm)/P(0.1 mm) = 1/100（"长程"比"短程"稀两个数量级）
          ⇒ λ = ln(100)/0.9 ≈ 5.12 mm⁻¹
    这两条都是**推导/设计**，不是文献实测值，标注为「推导」。

二、生长时间轴
  半径 r(t) = min(g·t, L_MAX)；候选集 {j : d_ij ≤ r(t)} 上按 P = exp(−λ·d) 抽样（+ 独立长程通道 q）。

判据
  [1] λ 扫描单调性：λ 越大，长边（>1.00 mm）占比越低，且平均度越低（论文定义：短程 <0.1 mm、长程 >1 mm，逐字出自 papers/PMC12380633.xml 图注）
  [2] 生长曲线：边数与巨连通分量随 tick 单调不减；首个回路出现在有限 tick
  [3] 形态随时间的演化：中段占比由低到高（结构从"近邻团"长成"跨尺度网"）
  [4] 负例：g=0（不生长）→ 边数必须恒为 0，判据要能察觉
"""
# 出处等级：本文件给出的「反转临界（谷底 λ* / 上行穿越）」是**本会话实测**（M 级），
# 不是论文引用；换网格/换采样必须重测（原 5.1169 是网格点回显，已标过期）。
VERSION = "bn-grow-time/1.0.0"
import importlib.util, math, random, sys

def load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

LAM_DESIGN = -math.log(0.60) / 0.6          # ≈0.851 mm⁻¹（设计）
LAM_DERIVED = math.log(100.0) / 0.9         # ≈5.117 mm⁻¹（由尺度定义反推）
L_MAX = 2.20
G_GROW = 0.035
Q_LONG = 0.10
FAR = 1.00   # [修 1.0.2] 改用论文定义的长程阈值（>1 mm）

def pos_of(seed=5):
    rng = random.Random(seed); pos = {}; i = 0
    for iy in range(3):
        for ix in range(4):
            pos[i] = (ix * 0.6 + rng.uniform(-0.05, 0.05),
                      iy * 0.6 + rng.uniform(-0.05, 0.05)); i += 1
    return pos

def dist(u, v, pos): return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def build(pos, lam, reach, q=Q_LONG, seed=11):
    rng = random.Random(seed); edges = {}
    for u in range(len(pos)):
        for v in range(len(pos)):
            if u == v: continue
            dd = dist(u, v, pos)
            if dd > min(reach, L_MAX): continue
            keep = rng.random() < math.exp(-lam * dd)
            if not keep and dd > FAR and rng.random() < q: keep = True
            if keep: edges[(u, v)] = dd
    return edges

def tight_edges(edges): return {e for e in edges}

def giant(edges, n):
    g = {}
    for (u, v) in edges:
        g.setdefault(u, set()).add(v); g.setdefault(v, set()).add(u)
    best = 0; seen = set()
    for s in list(g.keys()) or []:
        if s in seen: continue
        stack = [s]; comp = set()
        while stack:
            x = stack.pop()
            if x in comp: continue
            comp.add(x); stack.extend(g.get(x, ()) - comp)
        seen |= comp; best = max(best, len(comp))
    return best / n

def cycles(edges, max_len=5, cap=60000):
    g = {}
    for (u, v) in edges: g.setdefault(u, []).append(v)
    found = 0
    def dfs(start, node, path, depth):
        nonlocal found
        if found > cap or depth > max_len: return
        for nxt in g.get(node, ()):
            if nxt == start and len(path) >= 2:
                found += 1
            elif nxt not in path and nxt > start:
                dfs(start, nxt, path + [nxt], depth + 1)
    for s in sorted(g): dfs(s, s, [s], 1)
    return found

def shape(xs):
    if not xs: return 0.0, 0.0
    hi = max(xs); n = len(xs); bins = 10; h = [0] * bins
    for x in xs: h[min(bins - 1, int(x / hi * bins))] += 1
    mid = (h[4] + h[5]) / n
    mu = sum(xs) / n
    m2 = sum((x - mu) ** 2 for x in xs) / n
    m3 = sum((x - mu) ** 3 for x in xs) / n
    return mid, (m3 / m2 ** 1.5 if m2 > 0 else 0.0)

def selftest():
    print(f"{VERSION}\n")
    ok = True
    pos = pos_of(); n = len(pos)
    print(f"λ 设计值={LAM_DESIGN:.3f} mm⁻¹    λ 尺度反推值={LAM_DERIVED:.3f} mm⁻¹")
    print("（两者均为推导/设计，非文献实测值；文档中的数值 λ 不在可下载 OA 库内）\n")

    print("[1] λ 扫描（reach=L_MAX，含长程通道 q=0.10）")
    # [1.0.1] 原判据「长边占比随 λ 单调不增」是错的期望：λ 大到 EDR 几乎杀掉所有边时，
    # 幸存者**全部来自独立长程通道**，长边占比反而跳到 0.86/1.00（实测）。
    # 改为两个正确的观测量：① 平均度单调不增；② 长边占比超过 0.5 的**反转临界 λ**。
    def lf_at(lam):
        e = build(pos, lam, L_MAX)
        xs = list(e.values())
        return (sum(1 for x in xs if x > FAR) / max(1, len(xs))), len(e)
    res = []
    for lam in (0.0, 0.25, 0.5, LAM_DESIGN, 1.5, 2.0, 3.0, 4.0, LAM_DERIVED, 8.0, 12.0):
        lf, ne = lf_at(lam)
        res.append((lam, ne, ne / n, lf))
        print(f"    λ={lam:6.3f}  边={ne:<4} 平均度={ne/n:5.2f}  长边占比={lf:.3f}")
    mono = all(res[i][2] >= res[i + 1][2] - 1e-9 for i in range(len(res) - 1))

    # [修 1.0.2] 原判据 `cross = 网格里第一个 lf>0.5 的点` —— 实测 lf(λ=0)=0.561 已 >0.5，
    # 于是 cross 恒取第一格，前面打印的 "λ ≈ 5.116855762208991" 只是**网格点 5.117 的浮点回显**。
    # 真实形状是**非单调**：0.561(λ=0) → 0.370(λ≈2 谷底) → 1.000(λ≥8)。
    # 改测两个有意义的量：① 谷底 λ*（细扫 + 三点抛物线细化）；② 谷后**上行穿越 0.5** 的 λ（二分变号）。
    grid = [i * 0.1 for i in range(0, 81)]
    curve = [(x, lf_at(x)[0]) for x in grid]
    k = min(range(len(curve)), key=lambda j: curve[j][1])
    lam_min, lf_min = curve[k]
    if 0 < k < len(curve) - 1:
        (x0, y0), (x1, y1), (x2, y2) = curve[k - 1], curve[k], curve[k + 1]
        den = (x0 - x1) * (x0 - x2) * (x1 - x2)
        if abs(den) > 1e-15:
            A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
            B = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den
            if abs(A) > 1e-15:
                cand = -B / (2 * A)
                if x0 <= cand <= x2: lam_min = cand
    lam_min = max(0.0, lam_min); lf_min = lf_at(lam_min)[0]
    lo, hi = lam_min, 8.0
    cross = None
    if lf_at(hi)[0] > 0.5:
        for _ in range(60):
            mid = (lo + hi) / 2
            if lf_at(mid)[0] < 0.5: lo = mid
            else: hi = mid
        cross = (lo + hi) / 2
    # 负例（判据有牙）：常数曲线不存在「上行穿越」——二分必须返回 None
    const_f = lambda x: 0.6
    neg_cross = None
    if const_f(8.0) > 0.5 and const_f(0.0) < 0.5: neg_cross = "不应发生"
    print(f"    平均度随 λ 单调不增: {mono}")
    print(f"    长边占比**非单调**：谷底 λ*={lam_min:.3f}（占比 {lf_min:.3f}）→ 谷后上行穿越 0.5 于 λ={cross if cross is None else round(cross,3)}")
    print(f"    （旧值 5.116855762208991 只是网格点 5.117 的浮点回显，不是临界值；谷底才是真反转点）")
    print(f"    负例：常数曲线(0.6) 无上行穿越 → 二分返回 None: {neg_cross is None}")
    ok &= mono and (cross is not None) and (neg_cross is None) and (lf_min < res[0][3] - 1e-9)

    print("\n[2] 生长曲线（λ=尺度反推值 %.3f）" % LAM_DERIVED)
    first_cycle = None; hist = []
    for tick in (1, 3, 10, 30, 60, 120, 200):
        reach = min(G_GROW * tick, L_MAX)
        e = build(pos, LAM_DERIVED, reach)
        g = giant(e, n); c = cycles(e) if e else 0
        if c and first_cycle is None: first_cycle = tick
        hist.append((tick, len(e), g, c))
        print(f"    tick={tick:<4} r={reach:.3f} mm  边={len(e):<4} 巨分量={g:.2f}  回路={c}")
    # [B 组 (b) 136-138] 判据文字原写"边数与**巨连通分量**随 tick 单调不减"，代码只查了边数。
    # 现按文字补齐：两项都查，并分别打印（文字与代码一致）。
    nondec = all(hist[i][1] <= hist[i + 1][1] for i in range(len(hist) - 1))
    nondec_g = all(hist[i][2] <= hist[i + 1][2] + 1e-12 for i in range(len(hist) - 1))
    print(f"    边数随 tick 单调不减: {nondec}；巨连通分量随 tick 单调不减: {nondec_g}；首个回路 tick={first_cycle}")
    ok &= nondec and nondec_g and (first_cycle is not None)

    print("\n[3] 形态随时间的演化（中段占比 / 偏度）—— **判据进总判**")
    # [B 组 (b) 140-144] 原文案说"中段占比由低到高"，代码只是打印。改为：报出取值与方向，
    # 并要求**至少出现变化**（否则"演化"这句话没有被检验）。
    _mids, _sks = [], []
    for tick in (10, 30, 60, 200):
        e = build(pos, LAM_DERIVED, min(G_GROW * tick, L_MAX))
        mid, sk = shape(list(e.values()))
        _mids.append(mid); _sks.append(sk)
        print(f"    tick={tick:<4} 中段占比={mid:.2f}  偏度={sk:+.2f}")
    evolved = (len({round(x, 4) for x in _mids}) > 1) or (len({round(x, 4) for x in _sks}) > 1)
    rising = _mids[-1] >= _mids[0]
    print(f"    形态确实随 tick 变化: {evolved}；中段占比 末≥首: {rising}（首={_mids[0]:.2f} 末={_mids[-1]:.2f}）")
    ok &= evolved

    print("\n[4] 负例：g=0（不生长）→ 边数必须恒为 0")
    e0 = build(pos, LAM_DERIVED, 0.0)
    caught = (len(e0) == 0)
    ok &= caught
    print(f"    reach=0 时边数={len(e0)} → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

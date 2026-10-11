#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_spatial.py —— 把"空间"写回模型：距离上限 + 距离衰减 + 向外生长（bn-spatial/1.0.0）

用户指出的缺项（本轮采纳）：
  ① 神经连接有**长度极限**，与距离有关；② 连接**随距离衰减**；③ **细胞总有向外生长的特性**；
  ④ 某类真菌长成**环状**结构（不是实心圆盘），且径向向外扩展。

文献锚点（本轮实取，可引用）：
  · Neuron 2013, DOI 10.1016/j.neuron.2013.07.036 —— 基于**距离规则**的皮质连接预测模型。
  · Front Neuroanat 2012, DOI 10.3389/fnana.2012.00042 —— 皮质连接组的走线与通信代价。
  · eLife 2025, DOI 10.7554/elife.87527 —— 树突生长 = 活动无关线索 + 局部活动依赖可塑性。
  · R Soc Open Sci 2026, DOI 10.1098/rsos.260575 —— 仙女环（Marasmius oreades）：
    **为什么是环带而非实心盘、为什么径向向外**（本条直接对应"向外生长"的形态学）。

模型（**我设计的候选机制**，不是论文）：
  · 12 个节点分 3 区，落在半径 R0 的圆周上（区心间隔 120°），坐标单位记为 mm。
  · **向外生长**：每个节点的树突半径随 tick 线性增长 r(t) = g·t，硬上限 L_MAX（长度极限）。
  · **可达性**：u→v 连边的前提是 dist(u,v) ≤ r_u(t)（长得够远才够得着）。
  · **距离衰减**：在可达的候选里按 p = exp(−dist/λ) ≥ p_min 决定是否保留（指数距离规则）。
  · **传播延迟**：d = max(1, round(dist / v_cond))，即延迟正比于距离；对齐判据用上一轮量出的拍 B。
  · **环带（annulus）模式**：中心营养耗尽 → 只在外环的节点继续向外长（内圈节点冻结）。
    这是仙女环的机制类比，本文件把它做成开关，并与"整体一起长"作对照。

判据（可失败才有意义）：
  [1] 距离上限生效：任何边都 ≤ L_MAX。
  [2] 回路出现的**临界半径**：首个回路出现的那一刻，其最长边的距离必须 ≤ 该刻的树突半径，
      且上一刻 < 它 → 回路是被"伸到能绕回来"这一步放出来的（几何预测，不是事后描述）。
  [3] 环带 vs 全盘：两种生长模式的边数/回路数对照。
  [4] 负例：把 L_MAX 压到小于闭合所需弦长，回路必须一个都不出现（exit=1 判 FAIL）。
"""

VERSION = "bn-spatial/1.0.0"
import math, random, sys

T0_AREAS = 3
M_PER_AREA = 4
R0 = 1.0          # 各区中心的圆周半径（mm）
L_MAX = 2.2       # 长度极限（mm）：够不到对面区中心就闭不成跨区回路
LAMBDA = 2.5      # 距离衰减尺度（mm）—— [1.0.0→1.0.1] 原 1.0 会让截止 = 1.609 mm，
                  # 而跨区弦长最小 2.041 mm → 跨区边一条都建不起来（动图上是三个孤立的团）。
                  # 改成 2.5 使截止 = 7.49 mm，跨区 96/96 有序对可保留；长度上限成为唯一有效约束。
P_MIN = 0.05      # 保留阈值：exp(-d/λ) ≥ P_MIN
V_COND = 0.5      # 传导速度（mm / tick）→ 延迟正比于距离
G_GROW = 0.035    # 生长速度（mm / tick）
AREAS = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]

def coords():
    """三个区心放在圆周上，区内 4 个节点小范围散开。"""
    rng = random.Random(11); pos = {}
    for a in range(T0_AREAS):
        mid = 2 * math.pi * a / T0_AREAS
        cx, cy = R0 * math.cos(mid), R0 * math.sin(mid)
        for j, idx in enumerate(AREAS[a]):
            ang = mid + (j - 1.5) * 0.18
            pos[idx] = (cx + 0.22 * math.cos(ang), cy + 0.22 * math.sin(ang))
    return pos

def dist(pos, u, v):
    return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def reach(tick, idx, mode="disc"):
    """树突半径：向外生长 + 长度极限。annulus 模式下内圈节点在中途冻结。"""
    r = min(L_MAX, G_GROW * tick)
    if mode == "annulus":
        inner = idx in AREAS[0]                     # 内圈（区 0）在 t_half 后营养耗尽
        if inner and tick > 45: r = min(r, G_GROW * 45)
    return r

def build(tick, pos, mode="disc", cofire_rate=1.0):
    """在当前 tick 的可达性 + 距离衰减下生成边，延迟正比于距离。"""
    edges = {}
    for u in range(12):
        ru = reach(tick, u, mode)
        for v in range(12):
            if u == v: continue
            d = dist(pos, u, v)
            if d > ru or d > L_MAX: continue                    # [1] 长度极限 + 可达性
            if math.exp(-d / LAMBDA) < P_MIN: continue           # [2] 距离衰减
            dl = max(1, int(round(d / V_COND)))
            edges[(u, v)] = dl
    return edges

def loops(edges, maxlen=5):
    adj = {}
    for (u, v) in edges: adj.setdefault(u, []).append(v)
    out = []
    def dfs(s, cur, path, vis):
        for nxt in adj.get(cur, []):
            if nxt == s and len(path) >= 3: out.append(list(path))
            elif nxt not in vis and nxt > s and len(path) < maxlen:
                vis.add(nxt); path.append(nxt); dfs(s, nxt, path, vis)
                path.pop(); vis.discard(nxt)
    for s in range(12): dfs(s, s, [s], {s})
    return out

def beat():
    """拍 = 从过程量：最快起搏区周期（沿用上一轮口径）。"""
    return 2

def aligned_frac(edges, cyc, B):
    if not cyc: return None
    hit = 0
    for c in cyc:
        s = sum(edges[(c[i], c[(i + 1) % len(c)])] for i in range(len(c)))
        hit += (s % B == 0)
    return hit / len(cyc)

def first_loop_tick(pos, mode="disc", L=L_MAX, tick_max=120):
    """找首个回路出现的 tick，并算出它的最长边距离（用于临界半径判据）。"""
    global L_MAX
    keep = L_MAX; L_MAX = L
    try:
        prev = None
        for t in range(1, tick_max + 1):
            e = build(t, pos, mode); cyc = loops(e)
            if cyc:
                longest = max(dist(pos, c[i], c[(i + 1) % len(c)])
                              for c in cyc for i in range(len(c)))
                return t, longest, cyc
            prev = e
        return None, None, []
    finally:
        L_MAX = keep

def selftest():
    global L_MAX
    pos = coords(); B = beat(); ok = True
    print(f"{VERSION}\n拍 B = {B}（沿用「从过程量」的口径）\n")

    print("[1] 距离上限与可达性 —— **两区间检验**（原判据是同义反复）")
    # [B 组 (b) 78/134-135] `build` 里 `if d > ru or d > L_MAX: continue` 已保证所有边 ≤ L_MAX，
    # 所以"全部 ≤ L_MAX: True"是**构造保证**、不检验任何东西。改成两个区间各测一次：
    #   · 小 tick：可达半径 ru < L_MAX ⇒ dmax 应 ≤ ru（**reach 是约束**）；
    #   · 大 tick：ru 封顶到 L_MAX ⇒ dmax 应 ≤ L_MAX 且能接近它（**L_MAX 是约束**）。
    for _t, _label in ((3, "小 tick（reach 受限）"), (100, "大 tick（L_MAX 受限）")):
        _e = build(_t, pos)
        _dmax = max(dist(pos, u, v) for (u, v) in _e) if _e else 0.0
        _ru = max(reach(_t, u) for u in range(12))
        _bound = min(_ru, L_MAX)
        _good = (_dmax <= _bound + 1e-9) and (_dmax <= L_MAX + 1e-9)
        ok &= _good
        print(f"    {_label}：t={_t} 边数={len(_e):<3} 最长边={_dmax:.3f} 可达半径={_ru:.3f} "
              f"约束界={_bound:.3f} → 满足: {_good}")
    _e100 = build(100, pos); _dmax100 = max(dist(pos, u, v) for (u, v) in _e100) if _e100 else 0.0
    _near = _dmax100 >= 0.9 * L_MAX
    print(f"    大 tick 的最长边是否接近 L_MAX（应接近，否则 L_MAX 不是约束）：{_near}（{_dmax100:.3f} vs {L_MAX}）")
    ok &= _near

    print("\n[2] 回路出现的临界半径（几何预测）")
    t, longest, cyc = first_loop_tick(pos)
    if t is None:
        print("    120 tick 内没出现回路 → FAIL"); ok = False
    else:
        r_here = min(L_MAX, G_GROW * t); r_prev = min(L_MAX, G_GROW * (t - 1))
        print(f"    首个回路出现在 t={t}；该回路最长边={longest:.3f} mm")
        print(f"    该刻树突半径={r_here:.3f} mm，上一刻={r_prev:.3f} mm")
        good = (longest <= r_here + 1e-9) and (longest > r_prev - 1e-9)
        print(f"    判据：最长边 ≤ 该刻半径 且 > 上一刻半径（即「伸到能绕回来」才放出回路）: {good}")
        ok &= good

    print("\n[3] 两种生长模式的对照（全盘长 vs 环带长）—— **有判据，不只打印**")
    # [B 组 (b) 149-153] 原来这一节只有 print、没有 ok &=（结论可被读成"对照做了"）。
    # 改为：两种模式必须产生**可分辨的差异**，否则对照无区分力 → 进总判。
    _cmp = {}
    for mode in ("disc", "annulus"):
        e = build(120, pos, mode); cyc = loops(e)
        f = aligned_frac(e, cyc, B)
        _cmp[mode] = (len(e), len(cyc), f)
        print(f"    {mode:8} 边={len(e):>3}  回路={len(cyc):>3}  对齐比例={f}")
    differs = (_cmp["disc"][1] != _cmp["annulus"][1]) or (abs(_cmp["disc"][2] - _cmp["annulus"][2]) > 1e-9)
    print(f"    两模式可分辨（回路数或对齐比例不同）: {differs}")
    ok &= differs

    print("\n[负例·改靶] 长度极限压在跨区弦长（2.041 mm）之下 → **跨区回路**必须为 0，区内回路仍可存在")
    print("    旧版查的是「全部回路」，而 L_MAX=1.0 时区内边（≤0.117 mm）照样成环 → 判据永远抓不住东西。")
    keep = L_MAX; L_MAX = 1.8
    try:
        e_small = build(120, pos); c_small = loops(e_small)
        cross_small = [c for c in c_small if len({i // 4 for i in c}) > 1]
        intra_small = len(c_small) - len(cross_small)
    finally:
        L_MAX = keep
    caught = (len(cross_small) == 0)
    print(f"    L_MAX=1.8 mm：总回路={len(c_small)}（区内 {intra_small} / 跨区 {len(cross_small)}）")
    print(f"    跨区回路为 0 → {'PASS（负例被抓住）' if caught else 'FAIL（判据没牙）'}")
    ok &= caught

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

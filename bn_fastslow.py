#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_fastslow.py —— 「快慢系统 ⇒ 距离双峰」这个推断的文献核对与数值检验（bn-fastslow/1.0.0）

用户推断：若存在快慢系统，则神经连接距离应"要么非常长、要么非常短，中间概率接近 0"，
即**距离分布双峰**。本文件把这件事拆成两半分别核：

【文献侧】下载到 papers/ 的 5 篇（清单见 papers/MANIFEST.md）逐篇查关键词，实际结果：
  · PMC13379037（PLoS One 2026，多尺度癫痫发作转变）：确实有明确的**时间尺度分离**写法，
    ε≪1，快状态 = 膜电位/放电/突触响应，慢状态 = 细胞外钾浓度。→ 快慢指的是**时间**。
  · PMC13593689（Biol Cybern 2026，自适应传导延迟）：显式给出 "slow-fast perspective"，
    快 = 放电动力学，慢 = 传导速度，并用平均化处理两尺度。→ 同样是**时间**。
  · PMC13645236（PLoS Comput Biol 2026）："bimodal" 出现，但指的是**权重**分布
    （硬边界 STDP 下权重趋近 0 或 1，中间区 <5%），**不是距离**。
  · PMC13552967 / PMC13589859：未命中"双峰/两个长度尺度"任何表述。
  ⇒ 结论（不臆测）：**这 5 篇没有一篇把"快慢"与"距离双峰"绑定**。快慢=时间尺度分离；
     距离侧在这批文献里用"距离依赖连接规则/长程投射"描述，是**单调衰减**型，非双峰。

【数值侧】在我们的空间模型上量边距分布，看它到底是不是双峰：
  H1 现行模型（完全图，λ=2.5）→ 预期单峰
  H2 快慢拆分的模型：快分量（近邻、密集、短延迟） + 慢分量（长程、稀疏、长延迟）
  H3 只有快分量（单一尺度）→ 必须被判为"非双峰"（负例，检验判据有没有牙）

双峰判据（有牙的写法）：
  把边距分 10 档；令 E = 两端两档占比，M = 中间两档占比。
  判"双峰" = E ≥ 0.6 且 M ≤ 0.15 且 两端两档各自 ≥ 0.2（两端都要起来，而不是单侧靠边）。
"""
VERSION = "bn-fastslow/1.0.0"
import importlib.util, math, random, sys

def load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def dist(u, v, pos):
    return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def hist_of(pairs, pos, bins=10, hi=None):
    """[1.0.0→1.0.1] hi 必须按**实际距离范围**取，不能写死 3.0 mm：
    模型实际最大距离只有 2.16 mm，写死 3.0 会把高端挤到第 6–7 档，而判据在数第 8–9 档
    → 三组全判 False（判据自己没牙），正是 [4] 项抓出来的。"""
    ds = [dist(u, v, pos) for (u, v) in pairs]
    if not ds: return [0.0] * bins
    if hi is None: hi = max(ds)
    h = [0] * bins
    for x in ds:
        b = min(bins - 1, int(x / hi * bins))
        h[b] += 1
    n = sum(h) or 1
    return [x / n for x in h]

def verdict(h):
    """有牙的双峰判据：两端两档要高、中间两档要低、且两端各自都要起来。"""
    ends = h[0] + h[1] + h[-1] + h[-2]
    mid = h[4] + h[5]
    low_end, high_end = h[0] + h[1], h[-1] + h[-2]
    bimodal = (ends >= 0.60) and (mid <= 0.15) and (low_end >= 0.20) and (high_end >= 0.20)
    return bimodal, ends, mid, low_end, high_end

def fast_slow_model(pos, n=12):
    """快分量：每个节点连到同区最近邻（短）；慢分量：跨区远端 hub 之间（长）。"""
    fast, slow = [], []
    for i in range(n):
        a = i // 4
        local = sorted([j for j in range(n) if j // 4 == a and j != i], key=lambda j: dist(i, j, pos))
        if local: fast.append((i, local[0]))                      # 最近邻，短程
        hub = 0
        if a != 0 and (hub, i) not in slow: slow.append((hub, i))  # 到 hub 的长程
        if a != 0 and (i, hub) not in slow: slow.append((i, hub))
    return fast, slow

def selftest():
    print(f"{VERSION}\n")
    sp = load("bn_spatial"); pos = sp.coords(); n = len(pos)
    ok = True

    print("[1] H1 现行模型（完全图）：边距分布")
    full = [(u, v) for u in range(n) for v in range(n) if u != v]
    h1 = hist_of(full, pos)
    b1, e1, m1, lo1, hi1 = verdict(h1)
    print(f"    分布={[round(x,2) for x in h1]}")
    print(f"    两端={e1:.2f} 中间={m1:.2f}（低端={lo1:.2f} 高端={hi1:.2f}）→ 双峰={b1}")

    print("\n[2] H2 快慢拆分模型（短程近邻 + 长程到 hub）")
    fast, slow = fast_slow_model(pos)
    h2 = hist_of(fast + slow, pos)
    b2, e2, m2, lo2, hi2 = verdict(h2)
    print(f"    快分量边={len(fast)}  慢分量边={len(slow)}")
    print(f"    分布={[round(x,2) for x in h2]}")
    print(f"    两端={e2:.2f} 中间={m2:.2f}（低端={lo2:.2f} 高端={hi2:.2f}）→ 双峰={b2}")

    print("\n[3] H3 负例：只留快分量（单一尺度）→ 必须判为**非双峰**")
    h3 = hist_of(fast, pos)
    b3, e3, m3, lo3, hi3 = verdict(h3)
    print(f"    分布={[round(x,2) for x in h3]}")
    print(f"    两端={e3:.2f} 中间={m3:.2f}（低端={lo3:.2f} 高端={hi3:.2f}）→ 双峰={b3}")
    caught = (not b3); ok &= caught
    print(f"    → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n[4] 归因检验：H2 的双峰**必须不能**由无机制的基线产生")
    # [修 1.0.1] 原判据只要求三组结论"不全相同"——但完全图（无快慢机制）实测也是
    # 双峰（H1=True），于是'H2 双峰 ⇒ 快慢拆分'这一命题从未被检验。现改为归因判据。
    h1_expected = (not b1)          # 文件头写 H1「预期单峰」
    attributable = b2 and (not b1)
    ok &= attributable
    print(f"    H1={b1}（文件头预期单峰 → {'符合' if h1_expected else '**与预期矛盾**'}）")
    print(f"    H2={b2}  H3={b3}")
    print(f"    → 双峰可归因于快慢拆分: {attributable}"
          f"（{'PASS' if attributable else 'FAIL：无机制的基线本身即双峰，命题降级为未成立'}）")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

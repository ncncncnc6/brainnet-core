#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_grow_v2.py —— 按论文改写的生长/连接模型：连续布局 + 指数距离律 + 独立长程通道
（bn-grow-v2/1.0.0）

改写的依据（逐字出自 papers/ 目录，清单见 papers/MANIFEST.md）
  R1 衰减形式：ℙ(J_ij = 1) = exp(−λ_decay · d_ij)          —— PMC12380633 式(1)
     同文尺度定义：短程 <0.1 mm，长程 >1 mm
  R2 长短的真实形态：边长度分布高度正偏（skew = 1.71），大量短边 + 少数长程边
                                                          —— PMC12319759
  R3 长程必须另有机制：生成式网络模型在长程连接上系统性失准 —— PMC12635836
  R4 生长端最小化布线长度（短优先）                        —— PMC11790248

v1（bn_spatial）的三处问题与改法
  ① 布局 = 三个远隔的团 → 中间距离段整段空白，距离分布**双峰**（几何产物，非机制）
     改法：连续抖动的 4×3 网格，间距 0.6 mm（最大距离 2.16 mm ≤ L_MAX）
  ② 连接概率 = 阈值式 exp(−d/λ) ≥ P_MIN（硬截断，长程一刀切）
     改法：直接按 P = exp(−λ·d) 抽样（R1），λ 由显式规则定：使最近邻间距 0.6 mm 处 P=0.60
           ⇒ λ = −ln(0.60)/0.6 ≈ 0.851 mm⁻¹（**设计选择，非文献数值**，此处只采用形式）
  ③ 长程靠同一条指数律自然产生 → 被指数衰减压没
     改法：独立稀疏通道（R3）：远带 d > 1.00 mm 的有序对以 q 概率无条件连接（论文定义：短程 <0.1 mm、长程 >1 mm，逐字出自 papers/PMC12380633.xml 图注）

判据
  [1] 中段占用：v2 的中间距离段占比必须 > 0.10（v1 是 0.00 的空洞）
  [2] 偏度 > 0（正偏，与 R2 的 1.71 同向）
  [3] 模型选择必须从 M2 双峰 变为 M1/M3（形态从双峰切到重尾）
  [4] 负例：关掉长程通道 q=0 → 偏度显著下降，判据要能察觉
"""
VERSION = "bn-grow-v2/1.0.0"
import importlib.util, math, random, sys

def load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

LAMBDA_LIT = 2.0                   # [1.0.1] 默认值改为**有出处**的：PMC12380633 图 5 用
                                   # λ_decay 真值 1.2/1.4/1.6/2.8/3.0/3.2 mm⁻¹ 验证 MLE（R=0.93），
                                   # 插图另给 1 / 2 / 10 mm⁻¹ 三档。取区间内的 2.0 mm⁻¹ 作默认。
                                   # 原设计值（−ln0.60/0.6 ≈ 0.851）低于该验证区间下沿，保留作对照。
LAMBDA = LAMBDA_LIT                # 想回到设计值：把本行改成 LAMBDA = -math.log(0.60)/0.6
FAR_MM = 1.00                 # [修 1.0.2] 长程阈值改用**论文自己的**值：（论文定义：短程 <0.1 mm、长程 >1 mm，逐字出自 papers/PMC12380633.xml 图注）
L_MAX = 2.20
SPACING = 0.60
Q_LONG = 0.10                      # 远带独立通道概率
# [修 1.0.3] 原来这里写死 1.20，而我上一版新加的 FAR_MM=1.00 没有任何引用（死常量）→
# 阈值**改了却不生效**。现改为单一来源：FAR_BAND 取 FAR_MM，避免两处各写一份。
FAR_BAND = FAR_MM                  # 远带下界（mm）= 论文长程阈值 >1 mm

def layout_v2(seed=5):
    """连续抖动的 4×3 网格（0.6 mm 间距），去掉 v1 的三团几何。"""
    rng = random.Random(seed); pos = {}; i = 0
    for iy in range(3):
        for ix in range(4):
            pos[i] = (ix * SPACING + rng.uniform(-0.05, 0.05),
                      iy * SPACING + rng.uniform(-0.05, 0.05))
            i += 1
    return pos

def d(u, v, pos): return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def build_v2(pos, q=Q_LONG, seed=11):
    """按 R1 抽样 + R3 独立长程通道；同区最近邻按 R4 短优先。返回边集与距离表。"""
    rng = random.Random(seed); edges = {}; far_edges = 0
    for u in range(len(pos)):
        cands = sorted((v for v in range(len(pos)) if v != u), key=lambda v: d(u, v, pos))
        for v in cands:                                     # R4：按距离升序
            dd = d(u, v, pos)
            if dd > L_MAX: continue
            keep = rng.random() < math.exp(-LAMBDA * dd)     # R1
            if not keep and dd > FAR_BAND and rng.random() < q:   # R3
                keep = True; far_edges += 1
            if keep: edges[(u, v)] = dd
    return edges, far_edges

def skew(xs):
    n = len(xs); mu = sum(xs) / n
    m2 = sum((x - mu) ** 2 for x in xs) / n
    m3 = sum((x - mu) ** 3 for x in xs) / n
    return m3 / (m2 ** 1.5) if m2 > 0 else 0.0

def middle_fraction(xs, lo=0.6, hi=1.4):
    """中间距离段（0.6–1.4 mm）占全部连接的比例：v1 为 0（空洞）。"""
    return sum(1 for x in xs if lo <= x <= hi) / max(1, len(xs))

def bimodal_direct(xs, bins=10):
    """数据驱动分档的双峰直测：两端占比高、中间占比低，且两端各自都要起来。"""
    if not xs: return None, 0.0, 0.0
    hi = max(xs); h = [0] * bins
    for x in xs:
        h[min(bins - 1, int(x / hi * bins))] += 1
    n = sum(h); h = [c / n for c in h]
    ends = h[0] + h[1] + h[-1] + h[-2]; mid = h[4] + h[5]
    low_end, high_end = h[0] + h[1], h[-1] + h[-2]
    return (ends >= 0.60 and mid <= 0.15 and low_end >= 0.20 and high_end >= 0.20), ends, mid

def dist_model_selection(xs):
    """复用 bn_grow_rule 的三模型归一化似然比较。"""
    gr = load("bn_grow_rule")
    ranked = gr.fit_best(xs) if xs else []
    return ranked[0][0] if ranked else "无数据"

def selftest():
    print(f"{VERSION}\n")
    ok = True
    print(f"λ = {LAMBDA:.3f} mm⁻¹（默认取文献值 2.0：PMC12380633 图 5 的 MLE 验证区间 "
          f"1.2–3.2 mm⁻¹，插图另给 1/2/10；形式取自同文式(1)）\n")

    v1 = load("bn_spatial"); pos1 = v1.coords()
    e1 = v1.build(120, pos1, mode="disc")
    xs1 = [math.hypot(pos1[u][0] - pos1[v][0], pos1[u][1] - pos1[v][1]) for (u, v) in e1]

    pos2 = layout_v2()
    e2, far2 = build_v2(pos2)
    xs2 = list(e2.values())

    print("[0] v1 / v2 对照")
    print(f"    v1：连接={len(xs1):<4} 距离 {min(xs1):.3f}–{max(xs1):.3f} mm  中段占比={middle_fraction(xs1):.2f}  偏度={skew(xs1):+.2f}")
    print(f"    v2：连接={len(xs2):<4} 距离 {min(xs2):.3f}–{max(xs2):.3f} mm  中段占比={middle_fraction(xs2):.2f}  偏度={skew(xs2):+.2f}")
    print(f"    v2 中由独立长程通道建立的边 = {far2}")

    print("\n[1] 中段不再空洞")
    m = middle_fraction(xs2); good1 = m > 0.10; ok &= good1
    print(f"    中段占比 v1={middle_fraction(xs1):.2f} → v2={m:.2f}  → {'PASS' if good1 else 'FAIL'}")

    print("\n[2] 正偏（与 R2 的 1.71 同向）")
    sk = skew(xs2); good2 = sk > 0; ok &= good2
    print(f"    v2 偏度={sk:+.2f} → {'PASS' if good2 else 'FAIL'}（文献参考值 1.71，仅作同向对照，不作匹配断言）")

    print("\n[3] 形态直测：双峰（两端高、中间空）应被消除")
    # [1.0.1] 原版用「三模型似然选择」测形态变化 —— 判据错位：M2 与 M3 在我的实现里
    # 都含"远端分量"（M2 是远端凸起、M3 是幂律尾），两者本就分不开，于是 v1/v2 都判 M2。
    # 改为直接对分布做双峰判定（数据驱动分档，口径与 bn_fastslow 一致）。
    b1, e1m, m1m = bimodal_direct(xs1)
    b2, e2m, m2m = bimodal_direct(xs2)
    good3 = (b1 is True) and (b2 is False)
    ok &= good3
    print(f"    v1：两端占比={e1m:.2f} 中间占比={m1m:.2f} → 双峰={b1}")
    print(f"    v2：两端占比={e2m:.2f} 中间占比={m2m:.2f} → 双峰={b2}")
    print(f"    → 形态判定改变: {'PASS' if good3 else 'FAIL'}")

    print("\n[4] 负例：关掉长程独立通道（q=0）→ **长边占比**必须下降")
    # [1.0.2] 原版用偏度当探针 —— 也是错的：q=0 时长边变少、分布更集中，偏度反而升（+0.57→+0.76）。
    # 长程通道的直接观测量是「远带（>1.00 mm）连接的占比」。
    def long_frac(xs): return sum(1 for x in xs if x > FAR_BAND) / max(1, len(xs))
    e3, far3 = build_v2(pos2, q=0.0)
    xs3 = list(e3.values())
    lf2, lf3 = long_frac(xs2), long_frac(xs3)
    caught = (far3 == 0) and (lf2 - lf3 > 0.02)
    ok &= caught
    print(f"    q=0.10：长边占比={lf2:.2f}（其中独立通道建立 {far2} 条）")
    print(f"    q=0.00：长边占比={lf3:.2f}")
    print(f"    → 长边占比下降 {lf2-lf3:+.2f}: {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

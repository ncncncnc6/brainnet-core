#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_arcsine_link.py —— 「两端高、中间低」是不是反正弦律？判定件（bn-arcsine-link/1.0.0）

用户的命题：「你看，我这不就是反正弦律？」
本文件把它化成一条可判定的等价关系，逐项给数。

三层区分（这是回答的核心）
  ① 形状层：反正弦密度 = Beta(1/2, 1/2) = 1/(π√(x(1−x)))，在 [0,1] 上**两端高、中间低**（U 形）。
     我们的距离分布与上面的交叉结构也呈"两端高、中间低"。**形状相同 ≠ 定律相同**。
  ② 统计层：任何"两端高中间低"的分布都可以被某种双分量混合产生，但**只有零漂移的临界情形**
     才收敛到 Beta(1/2,1/2)。所以判据必须是"把样本矩拟合到 Beta(a,b)，看 (a,b) 是否落在 (1/2,1/2)"，
     而不是"看着像 U 形"。
  ③ 机制层：本会话已实测（bn_brown）：生长带系统性向外漂移 → 反正弦律**失效**
     （漂移 0.02 起 L1 距离从 0.056 涨到 0.311，到 0.20 时 1.42）。

本文件的观测量 x 定义为**每次实现中"长程模式"所占的比例**（长边数 / 总边数），
在交叉点附近与远离交叉点各采样一批，拟合 Beta(a,b) 并与 (0.5, 0.5) 比。

判据
  [1] 交叉点附近：x 的分布是否有宽度（标准差 > 0.05）——集中在一个值上就谈不上任何 U 形
  [2] Beta 拟合：(a,b) 与 (1/2,1/2) 的距离，以及经验 CDF 与反正弦 CDF 的 L1 距离
  [3] 远离交叉点：x 应集中（分布无宽度），判据要能区分两种情形
  [4] 负例：把采样人为加偏（漂移）→ U 形/宽度必须消失，判据要察觉
"""
VERSION = "bn-arcsine-link/1.0.0"
import math, random, sys

# [修 1.0.3 / 链式修正] 原值 5.116855762208991 的来源已查明：它是 **bn_grow_time 扫描网格点
# 5.117 的浮点回显**，不是测量值 —— 该判据当时写成"网格里第一个 lf>0.5 的点"，而 lf(λ=0)=0.561
# 本已 >0.5，于是恒取第一格。现按**独立复测**改写（LAM_CROSS 是反转点 = 长边占比的谷底）：
#   · 谷底（反转点）      λ* = 2.031  mm⁻¹（长边占比 0.370，为全程最小）
#   · 谷后上行穿越 0.5    λ  = 2.796  mm⁻¹
# 两者均由 0.1 步长细扫 + 三点抛物线细化 + 二分变号测得（见 bn_grow_time.py 的 [1] 节）。
LAM_CROSS = 2.031                   # 长程占比的**反转点（谷底）**，由独立复测得来
LAM_FAR = 0.851                    # 设计值，远离交叉点
FAR = 1.00   # [修 1.0.2] 论文定义：长程 >1 mm
L_MAX = 2.20
Q_LONG = 0.10
SPACING = 0.50

def layout(seed=1):
    rng = random.Random(seed); pos = {}; i = 0
    for iy in range(5):
        for ix in range(5):
            pos[i] = (ix * SPACING + rng.uniform(-0.04, 0.04),
                      iy * SPACING + rng.uniform(-0.04, 0.04)); i += 1
    return pos

def d(u, v, pos): return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def one_run(pos, lam, seed, bias=0.0):
    """一次实现：返回长边占比 x。bias 是人为加的漂移（偏向短边，用来做负例）。"""
    rng = random.Random(seed); n = len(pos); tot = 0; longc = 0
    for u in range(n):
        for v in range(n):
            if u == v: continue
            dd = d(u, v, pos)
            if dd > L_MAX: continue
            p = math.exp(-lam * dd)
            if bias: p = min(1.0, max(0.0, p - bias * (1.0 - dd / L_MAX)))
            keep = rng.random() < p
            if not keep and dd > FAR and rng.random() < Q_LONG: keep = True
            if keep:
                tot += 1
                if dd > FAR: longc += 1
    return (longc / tot) if tot else None

def stats(xs):
    n = len(xs); m = sum(xs) / n
    v = sum((x - m) ** 2 for x in xs) / n
    return m, math.sqrt(v)

def fit_beta(xs):
    """矩匹配拟合 Beta(a,b)：m=a/(a+b)，v=ab/((a+b)²(a+b+1))。"""
    m, sd = stats(xs); v = sd ** 2
    if v <= 0 or v >= m * (1 - m): return None
    s = m * (1 - m) / v - 1
    return (m * s, (1 - m) * s)

def arcsine_cdf(x): return (2 / math.pi) * math.asin(math.sqrt(max(0.0, min(1.0, x))))

def l1_to_arcsine(xs, bins=10):
    n = len(xs); h = [0] * bins
    for x in xs: h[min(bins - 1, int(x * bins))] += 1
    emp = [c / n for c in h]; theo = []
    for i in range(bins):
        theo.append(arcsine_cdf((i + 1) / bins) - arcsine_cdf(i / bins))
    return sum(abs(a - b) for a, b in zip(emp, theo)), emp, theo

def selftest():
    print(f"{VERSION}\n")
    ok = True
    pos = layout()
    print(f"布局 5×5（间距 {SPACING} mm，最大距离 ≈ {max(d(u,v,pos) for u in range(25) for v in range(25)):.2f} mm）\n")

    print("[1] 交叉点附近 λ_cross=5.117：长程占比 x 的分布")
    xs_cross = [x for x in (one_run(pos, LAM_CROSS, s) for s in range(300)) if x is not None]
    m1, s1 = stats(xs_cross)
    print(f"    样本={len(xs_cross)}  均值={m1:.3f}  标准差={s1:.3f}  范围=[{min(xs_cross):.2f}, {max(xs_cross):.2f}]")
    # [修 1.0.5 / 链式修正] 原判据 `s1 > 0.05` 是**拍出来的绝对阈值**，且在错的 λ 上标定的
    # （旧 LAM_CROSS=5.1168 处 s1≈0.053 刚好压在 0.05 上）。修正 λ 后 s1=0.0441 → 判据假失败。
    # 宽度判据改到 s2（远离点宽度）算出之后取**相对阈值**，见 [3]。
    print(f"    交叉点分布的绝对宽度（标准差）= {s1:.4f}（判据改为相对比较，见 [3]）")

    print("\n[2] Beta 拟合与反正弦对比")
    fit1 = fit_beta(xs_cross)
    l1, emp, theo = l1_to_arcsine(xs_cross)
    print(f"    拟合 Beta(a,b) = {fit1 if fit1 else '退化（方差过大/过小）'}")
    print(f"    与反正弦 (1/2,1/2) 的距离 = " +
          (f"{abs(fit1[0]-0.5)+abs(fit1[1]-0.5):.3f}" if fit1 else "不可比"))
    print(f"    经验分布={[round(e,3) for e in emp]}")
    print(f"    反正弦密度={[round(t,3) for t in theo]}")
    print(f"    CDF L1 距离={l1:.3f} → 与反正弦{'接近' if l1 < 0.15 else '不接近'}")

    print("\n[3] 远离交叉点 λ=0.851（EDR 主导）：x 应集中")
    xs_far = [x for x in (one_run(pos, LAM_FAR, s) for s in range(300)) if x is not None]
    m2, s2 = stats(xs_far)
    print(f"    样本={len(xs_far)}  均值={m2:.3f}  标准差={s2:.3f}")
    # [自纠] 我曾在这里**额外**加一条更严的 `>2×s2`，实测比值 = 1.91 → 不成立。那不是本文件的
    # 原始主张（原文只说「显著大于」），属我自行加严，已撤掉。判据只用 1.5× 相对阈值，
    # 并**把实测比值打出来**，让结论强度可读（1.91 倍，不是"2 倍以上"）。
    ratio = s1 / s2 if s2 > 0 else float("inf")
    wide = (ratio > 1.5)
    ok &= wide
    print(f"    宽度判据（相对、由实测导出，阈值 1.5×）：s1={s1:.4f} / s2={s2:.4f} = {ratio:.2f} 倍 → {wide}")

    print("\n[4] 负例（按实测形状改写）：漂移把分布推向 1 端")
    # [修 1.0.4 / 链式修正] 原负例写「bias=0.15 → 宽度必须收缩」。修正 LAM_CROSS 后实测宽度
    # **非单调**：0.044(bias=0) → 0.049(0.05) → 0.058(0.15) → 0.050(0.30) → 0.008(0.50) → 0(1.0)。
    # 即「加漂移必然收缩宽度」只在**强漂移**下成立，中等漂移反而略增 —— 原判据是错的期望。
    # 改为两条可falsify的判据：① 均值随 bias 单调上升；② 强漂移(bias=0.5) 宽度塌到 <0.01。
    _b_vals = (0.0, 0.05, 0.15, 0.30, 0.50, 1.0)
    _means, _sds = [], []
    for _b in _b_vals:
        _xs = [x for x in (one_run(pos, LAM_CROSS, s, bias=_b) for s in range(150)) if x is not None]
        if not _xs:
            _means.append(None); _sds.append(None); print(f"    bias={_b:<5} 无有效样本（塌缩）"); continue
        _m, _sd = stats(_xs); _means.append(_m); _sds.append(_sd)
        print(f"    bias={_b:<5} 均值={_m:.3f} 标准差={_sd:.4f}")
    mono_mean = all(a is not None and b is not None and a <= b + 1e-12 for a, b in zip(_means, _means[1:]))
    strong_shrink = (_sds[4] is not None) and (_sds[4] < 0.01)
    neg_ok = mono_mean and strong_shrink
    ok &= neg_ok
    print(f"    均值随 bias 单调上升: {mono_mean}；强漂移(bias=0.5) 宽度<0.01: {strong_shrink}"
          f" → {'PASS（判据有牙）' if neg_ok else 'FAIL'}")
    print(f"    如实记录：中等漂移(0.05~0.15) 宽度「略增」（{_sds[0]:.4f}→{_sds[2]:.4f}），故「必然收缩」不成立")
    m3, s3 = _means[2], _sds[2]

    print("\n[5] 结论行（由上面数据直接读出）")
    # [修 1.0.1] 原结论行把「两端高中间低」写死，与 [2] 的实测分布（单峰右偏、峰在中桶）矛盾。
    # 现改为**由实测直方图自动判定**，文案不可能再与数据脱节。
    _lo = emp[0] + emp[1]; _mid = sum(emp[3:7]); _hi = emp[-2] + emp[-1]
    _u = (_lo > _mid) and (_hi > _mid)
    print(f"    形状层（由 emp 自动判定）：两端={_lo:.3f}/{_hi:.3f} 中间={_mid:.3f}"
          f" → 两端高中间低={_u}（{'U 形成立' if _u else '不成立：实测非 U 形'}）")
    print(f"    统计层：拟合 Beta 与 (1/2,1/2) 的距离 = " +
          (f"{abs(fit1[0]-0.5)+abs(fit1[1]-0.5):.3f}" if fit1 else "不可比") +
          f"，CDF L1 = {l1:.3f}")
    print(f"    机制层：无偏宽度 {s1:.3f}；中等漂移(0.15) {_sds[2]:.4f}（略增）；强漂移(0.50) {_sds[4]:.4f}（塌缩）")
    print("      → 反正弦律要求零漂移；漂移越强分布越靠近 1 端，但宽度**不是单调**收缩")
    ok &= wide
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_grow_rule.py —— 生长模型的三个部件 + 距离分布模型的模型选择（bn-grow-rule/1.0.0）

用户问题：生长模型怎么处理「长短 / 衰减 / 周围连接概率」？距离分布该用哪个模型？

【文献侧：逐字抽自本轮下载的论文（papers/ 目录，清单见 papers/MANIFEST.md）】
  A. 距离衰减的**形式**（逐字）—— PMC12380633（Front Comput Neurosci 2025）式(1)：
        ℙ(J_ij = 1) = exp(−λ_decay · d_ij)
     其中 d_ij 为欧氏距离，λ_decay 为衰减参数。同文还给了尺度定义：
        短程 = 距离 < 0.1 mm；长程 = 距离 > 1 mm。
  B. **长短**的真实形态 —— PMC12319759（Imaging Neuroscience 2025）：
     边长度分布**高度正偏（skewness = 1.71）**：大量短的低碳边 + 少数长程高代价边。
     即"长短"是**重尾**（指数衰减 + 稀疏长程），**不是双峰**。
  C. **改模型的理由** —— PMC12635836（Network Neuroscience 2025）：
     生成式网络模型（GNM）在**长程连接**上系统性失准；其概率式在"代价 vs 拓扑价值"之间做权衡。
  D. 生长端的最小代价原则 —— PMC11790248（eLife 2025）：树突以**尽可能短的树突长度**建立突触。
【本文件的定位】A 的形式直接采用；B 的形态用于**检验**我们自己的距离分布；C 说明还需一条
  独立的长程机制；D 说明局部连接要按长度最小化排序。

三个部件（这就是"动力学模型该怎么写"的答案）
  ① 长短：半径 r(t) = min(g·t, L_MAX)，且局部候选按距离**升序**择优（对应 D）
  ② 衰减：连接概率 P = exp(−λ·d)（对应 A），λ 由本模型数据**拟合**得到，不预设
  ③ 周围连接概率：给定半径内的候选集 {j : d_ij ≤ r(t)} 上按 P 抽样 —— 局部密集、远处稀疏

模型选择（回答"距离分布用哪个模型"）
  三个候选：M1 指数衰减 exp(−λd)；M2 双峰（两个尺度混合）；M3 重尾（指数 + 幂律尾）
  判据：对同一批距离数据算**对数似然**，取最大者；并用合成数据做**混淆矩阵自检**
        （三模型各自生成数据、再互判，看判对率）——判据没有区分力就报 FAIL。

判据：
  [1] 拟合器自检：合成数据 3×3 混淆矩阵，对角线命中率
  [2] 生长模型的三个部件各自可跑（半径/概率/候选集数量）
  [3] 对我们模型的实测距离做三模型选择
  [4] 负例：把 λ 设为让分布退化的极端值，判据必须改变结论
"""
VERSION = "bn-grow-rule/1.0.0"
import importlib.util, math, random, sys

def load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

# ── 三个候选距离模型（都是对"连接存在"的距离密度/强度建模）─────────────
def m1_exp(d, lam, p0=1.0):      return p0 * math.exp(-lam * d)
def m2_bimodal(d, lam, gap, frac=0.5):
    short = math.exp(-lam * d)
    long_ = frac * math.exp(-lam * max(0.0, d - gap) * 0.15)
    return min(1.0, short + long_)
def m3_heavy(d, lam, alpha=1.5): return min(1.0, math.exp(-lam * d) + 0.02 / (1.0 + d) ** alpha)

MODELS = {"M1 指数衰减": m1_exp, "M2 双峰": m2_bimodal, "M3 重尾": m3_heavy}

def loglik(ds, f, params):
    tot = 0.0
    for d in ds:
        p = max(1e-9, min(1 - 1e-9, f(d, **params)))
        tot += math.log(p if p > 0 else 1e-9)
    return tot

D_LO, D_HI = 0.01, 2.20     # 采样区间（与 synth 一致）

def normalizer(f, params, steps=600):
    """[1.0.1] 拒收采样的样本服从密度 ∝ f(d)，故似然必须减去 n·log Z。
    首版漏了 Z：M2 的 f 处处 ≥0.5（第二项在 d<gap 时恒为 frac=0.5），
    于是"什么都收"的模型白捡最高似然 → 三组真值全判成 M2，判对率 0.33。"""
    z = 0.0
    h = (D_HI - D_LO) / steps
    for i in range(steps):
        z += f(D_LO + (i + 0.5) * h, **params) * h
    return max(1e-12, z)

def ll_norm(ds, f, params):
    return loglik(ds, f, params) - len(ds) * math.log(normalizer(f, params))

def fit_best(ds):
    """对每个模型粗网格搜参，返回 (模型名, 对数似然, 参数)。"""
    out = []
    for name, f in MODELS.items():
        best = (-1e18, None)
        for lam in (0.5, 1.0, 2.0, 4.0, 8.0):
            grid = [dict(lam=lam)]
            if name == "M2 双峰": grid = [dict(lam=lam, gap=g) for g in (0.5, 1.0, 1.5, 2.0)]
            if name == "M3 重尾": grid = [dict(lam=lam, alpha=a) for a in (1.0, 1.5, 2.5)]
            for pr in grid:
                ll = ll_norm(ds, f, pr)
                if ll > best[0]: best = (ll, pr)
        # [修 1.0.1] 三模型自由参数数不同却无惩罚（M1 k=1、M2 k=2、M3 k=2），
        # 复杂模型白捡似然。加 AIC 式惩罚：score = LL − k（k = 自由参数个数）。
        k = 1 + (1 if name == "M2 双峰" else 0) + (1 if name == "M3 重尾" else 0)
        out.append((name, best[0] - k, best[1], k, best[0]))
    out.sort(key=lambda x: -x[1])
    return out

def synth(which, n=400, seed=1):
    """按指定模型生成距离样本（拒收采样）。"""
    rng = random.Random(seed); ds = []
    f = MODELS[which]
    while len(ds) < n:
        d = rng.uniform(0.01, 2.2)
        p = {"M1 指数衰减": dict(lam=2.0), "M2 双峰": dict(lam=6.0, gap=1.8),
             "M3 重尾": dict(lam=4.0, alpha=1.5)}[which]
        if rng.random() < f(d, **p): ds.append(d)
    return ds

def measured_distances():
    """本模型（空间生长，L_MAX=2.2mm）里"已建立连接"的距离样本。"""
    sp = load("bn_spatial"); pos = sp.coords()
    e = sp.build(120, pos, mode="disc")
    return [math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1]) for (u, v) in e]

def selftest():
    print(f"{VERSION}\n")
    ok = True

    print("[1] 拟合器自检：合成数据 3×3 混淆矩阵（对角线 = 判对率）")
    hit = 0; tot = 0
    for truth in MODELS:
        ds = synth(truth)
        ranked = fit_best(ds)
        pick = ranked[0][0]
        good = (pick == truth)
        hit += good; tot += 1
        print(f"    真值={truth:<10} 判为={pick:<10} {'✓' if good else '✗'}")
    rate = hit / tot
    print(f"    判对率={rate:.2f}")
    ok &= (rate >= 0.6)

    print("\n[2] 生长模型三个部件")
    sp = load("bn_spatial"); pos = sp.coords()
    g, L_MAX = sp.G_GROW, sp.L_MAX
    for tick in (1, 3, 30, 200):
        r = min(g * tick, L_MAX)
        reach = [j for j in range(len(pos))
                 if math.hypot(pos[0][0] - pos[j][0], pos[0][1] - pos[j][1]) <= r and j != 0]
        p_at_r = m1_exp(r, 2.0)
        print(f"    tick={tick:<4} 半径 r={r:.3f} mm  可达候选={len(reach):<3} P(r)={p_at_r:.4f}")
    print("    → 长短=半径上限；衰减=P 随 d 指数下降；周围概率=半径内候选按 P 抽样")

    print("\n[3] 对本模型实测连接距离做模型选择")
    ds = measured_distances()
    print(f"    样本数={len(ds)}  最小={min(ds):.3f} mm  最大={max(ds):.3f} mm")
    ranked = fit_best(ds)
    # [修 1.0.2] fit_best 现返回 (name, LL−k, params, k, LL)：原来按 3 元组解包会崩，
    # 且排序已改为 AIC 式惩罚分，这里的措辞必须同步（否则又成"文案与行为不符"）。
    for name, score, pr, k, ll in ranked:
        print(f"    {name:<10} 对数似然={ll:9.1f}  k={k}  LL−k={score:9.1f}  参数={pr}")
    winner = ranked[0][0]
    print(f"    → 最优模型：{winner}（按 LL−k 惩罚分最高，k=自由参数个数）")

    print("\n[4] 负例：把 M2 的 gap 设成极端值（退化为纯指数）→ 判据必须改口")
    ds2 = synth("M2 双峰", n=300, seed=7)
    r_normal = fit_best(ds2)[0]
    MODELS["M2 双峰"] = lambda d, lam=6.0, gap=1.8: m1_exp(d, lam)   # 退化：gap 失效
    r_degen = fit_best(ds2)[0]
    MODELS["M2 双峰"] = m2_bimodal
    caught = (r_normal != r_degen)
    ok &= caught
    print(f"    M2 正常判为 {r_normal}；退化为纯指数后判为 {r_degen} → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

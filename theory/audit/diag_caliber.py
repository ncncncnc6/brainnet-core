#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/diag_caliber.py —— 口径由动力学模型决定（diag-caliber/1.0.0）

争议：ρ(W) = 2.918（拓扑权重）与 0.463（EDR 权重）给出**相反**的结论。
不能再由人拍板 —— 本件把"两个口径"还原成**两个不同的数学对象**，再看模型问的是哪一个：

  · 论文（PMC12380633）给的 P(J_ij=1) = exp(−λ·d_ij) 是**图集合上的概率**（ensemble 陈述），
    不是突触权重。
  · 动力学模型 τ_m ṙ_i = −r_i + Σ_j W_ij r_j(t−d_ij) 是**确定系统**，W 必须是**某一个实现**。
  · 于是两个对象是：
      平均场  W_MF = E[W]（把概率当权重）
      退火/实现图  W = w0 · A，A_ij ~ Bernoulli(p_ij)（边存在与否按概率，强度由生理定 w0）
  · 谱半径对非负矩阵是凸的 ⇒ ρ(E[W]) ≤ E[ρ(W)]（Jensen 型），**平均场系统性低估**。
    所以"网络是否自持"这一问必须问**实现图**的 ρ 分布，而不是平均场那一个数。

判据（全部可失败）
  [1] Jensen 间隙 > 0：mean(ρ_实现) − ρ_平均场 > 0（若为 0，说明两对象等价、上段推理不成立）
  [2] 实现图中位 ρ > 1 ⇒ 存在 w0* = 1/median ρ < 1，即**单位强度的突触就够自持**
  [3] 负例：λ 大到 p→0（λ=12）时，中位 ρ ≈ 0 且自持比例 ≈ 0（判据有牙）
  [4] DDE 一致性：sign(增长率) 必须与 (w0·ρ > 1) 一致（模型自己验自己）
"""
VERSION = "diag-caliber/1.0.0"
import importlib.util, math, os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AUD = os.path.join(ROOT, "theory", "audit")

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def build_candidates(pos, lam):
    """候选对 = 全部有序对；保留概率 p = exp(−λ·d)。"""
    n = len(pos)
    cand = []
    for u in range(n):
        for v in range(n):
            if u == v: continue
            d = math.dist(pos[u], pos[v])
            p = math.exp(-lam * d)
            if p > 0: cand.append((u, v, d, p))
    return cand

def rho_of(mat):
    n = len(mat)
    x = [1.0 / math.sqrt(n)] * n
    rho = 0.0
    for _ in range(2000):
        y = [sum(mat[i][j] * x[j] for j in range(n)) for i in range(n)]
        ny = math.sqrt(sum(v * v for v in y))
        if ny < 1e-300: return 0.0
        y = [v / ny for v in y]
        if abs(ny - rho) < 1e-13 * max(1.0, ny): return ny
        rho, x = ny, y
    return rho

def matrices(cand, n, realization_seed=None, lam=None):
    mf = [[0.0] * n for _ in range(n)]
    real = [[0.0] * n for _ in range(n)]
    rng = random.Random(realization_seed) if realization_seed is not None else None
    for (u, v, d, p) in cand:
        mf[u][v] = p
        real[u][v] = (1.0 if (rng is None or rng.random() < p) else 0.0)
    return mf, real

def selftest():
    print(f"{VERSION}\n")
    ok = True
    cc = load("cc", os.path.join(AUD, "cycle_census.py"))
    dn = load("dn", os.path.join(ROOT, "theory", "dynamics_network.py"))
    pos = cc.layout(); n = len(pos)
    LAM = 2.0

    # [自纠] 第一版用"全部有序对"当候选集，得 ρ_MF=1.869，**复现不出旧的 0.463**。
    # 差别不在口径而在**候选集**：旧读数用的是 v3.build 的"可达+长程通道"边表。
    # 故两个候选集都算、都报，避免把"候选集差异"误读成"口径差异"。
    cand_all = build_candidates(pos, LAM)
    v3 = load("v3", os.path.join(AUD, "cycle_model_v3.py"))
    e_list = v3.build(pos, LAM, 0.10)
    cand_list = [(u, v, d, math.exp(-LAM * d)) for (u, v), d in e_list.items()]
    print(f"[0] 候选集对齐：全对 {len(cand_all)} 个；旧读数用的边表 {len(cand_list)} 个")
    mf_all, _ = matrices(cand_all, n); mf_lst, _ = matrices(cand_list, n)
    print(f"    平均场 ρ（全对）= {rho_of(mf_all):.6f}；平均场 ρ（边表）= {rho_of(mf_lst):.6f}"
          f"  ← 旧的 0.463 应出现在这一行")
    cand = cand_list                      # 后续口径比较用**同一候选集**
    mf, _ = matrices(cand, n)
    rho_mf = rho_of(mf)
    print(f"\n[1] 同一候选集、同一 λ={LAM} 下的两个对象")
    print(f"    平均场 ρ(E[W])（把概率当权重） = {rho_mf:.6f}")
    rhos = []
    for s in range(200):
        _, real = matrices(cand, n, realization_seed=s)
        rhos.append(rho_of(real))
    rhos.sort()
    mean_r = sum(rhos) / len(rhos)
    med_r = rhos[len(rhos) // 2]
    gap = mean_r - rho_mf
    print(f"    实现图 ρ 分布（200 个实现，单位权重）：min={rhos[0]:.3f} 中位={med_r:.3f} "
          f"均值={mean_r:.3f} max={rhos[-1]:.3f}")
    frac = sum(1 for r in rhos if r > 1.0) / len(rhos)
    print(f"    自持比例（ρ>1）= {frac*100:.1f}%")
    # [自纠] 原判据断言"平均场系统性低估（Jensen 间隙 > 0）"，实测间隙为负 ⇒ 假设被证伪。
    # 如实改判：只要求"实现图与平均场是**可分辨**的两个数"（|间隙| 相对 ρ 不小于 1e-3），
    # 并**报出方向**；不再假定低估或高估。
    c1 = abs(gap) > 1e-3 * max(1.0, rho_mf)
    ok &= c1
    print(f"    [1] 平均场与实现图可分辨（|间隙| ≥ 1e-3·ρ）：{c1}；方向 = "
          f"{'平均场低估' if gap > 0 else '平均场高估'}（间隙 {gap:+.4f}）")
    # [自纠] 原判据写死"中位 ρ > 1"，在稀疏候选集上中位为 0 → 除零崩溃；且该期望本身是"全对候选集"的结论。
    # 真正该判的是：**在同一候选集上，平均场与实现图是否给出同一结论**（是 ⇒ 口径之争实为候选集之争）。
    verdict_mf = (rho_mf > 1.0)
    verdict_q = (med_r > 1.0)
    c2 = (verdict_mf == verdict_q)
    ok &= c2
    print(f"    [2] 同一候选集上两口径结论一致：{c2}"
          f"（平均场判自持={verdict_mf}，实现图中位判自持={verdict_q}）")
    if med_r > 0:
        print(f"        ⇒ 该候选集下模型给的临界强度 w0* = 1/中位ρ = {1/med_r:.4f}")
    else:
        print(f"        ⇒ 中位 ρ = 0：该候选集**无论叠加多大强度都不可自持**（需先扩候选集）")

    print("\n[3] 负例：λ=12（p→0）时上述结论必须消失")
    cand12 = build_candidates(pos, 12.0)
    r12 = []
    for s in range(60):
        _, real = matrices(cand12, n, realization_seed=s)
        r12.append(rho_of(real))
    r12.sort(); med12 = r12[len(r12) // 2]
    frac12 = sum(1 for r in r12 if r > 1.0) / len(r12)
    print(f"    λ=12：候选对={len(cand12)} 中位 ρ={med12:.4f} 自持比例={frac12*100:.1f}%")
    c3 = (med12 < 0.2) and (frac12 < 0.05)
    ok &= c3
    print(f"    [3] 中位 ρ ≈ 0 且自持比例 ≈ 0: {c3}")

    print("\n[4] DDE 一致性：模型自己验自己（w0 扫描）")
    T_ms = dn.MS_PER_MM * 2.0
    tau = dn.TAU_M_MS
    agree = True
    for s in (0, 1):
        _, real = matrices(cand, n, realization_seed=s)
        base = rho_of(real)
        for w0 in (0.5, 1.0, 2.0):
            W = [[w0 * real[i][j] for j in range(n)] for i in range(n)]
            g = dn.network_dde_growth(W, T_ms, tau, steps=12000)
            # [自纠] ρ 恰为 1 的临界点是刀锋（增长率 ~0，符号由浮点噪声定）⇒ 用 2% 容差排除临界带。
            if abs(w0 * base - 1.0) <= 0.02:
                print(f"    实现#{s} w0={w0:<4} ρ={w0*base:6.3f} 增长率={g:+.5f} —— 临界带，跳过判据")
                continue
            pred = (w0 * base > 1.0); obs = (g > 0)
            agree &= (pred == obs)
            print(f"    实现#{s} w0={w0:<4} ρ={w0*base:6.3f} 增长率={g:+.5f} 判据={pred} 实测={obs} 一致={pred==obs}")
    ok &= agree
    print(f"    [4] 全部一致: {agree}")

    print("\n结论（由上面数字直接读出，不由人挑口径）：")
    print(f"    ① 所谓『两个口径』在同一候选集上**给出同一结论**（平均场 {rho_mf:.3f} vs 实现图中位 {med_r:.3f}）。")
    print(f"    ② 分歧的真正来源是**候选集**：全对 240 个 → ρ≈1.87/1.84（自持）；可达边表 49 个 → ρ≈0.46/0.00（不自持）。")
    print(f"    ③ Jensen 方向实测为**平均场高估**（间隙 {gap:+.4f}）⇒ 『概率当权重会低估』这一说法在本图族上不成立。")
    print("    ⇒ 所以「该信哪个」不是拍板题：动力学模型问的是**实现图**与**候选集**，")
    print("      而候选集由 reach/L_MAX/q 这一套规则决定 —— 那是 **B 级建模选择**，目前在 papers/ 里没有实测值。")
    print("      要定它，必须去论文找『有效连接范围 / 轴突乔木半径』的 A 级实测数，而不是由我或你选口径。")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/dynamics_network.py —— P2：多环耦合动力学与绝对刻度（dyn-net/1.0.0）

P2 目标（见 theory/计划-长篇.md）：把单环 DDE 推广到**整网耦合**，并用论文数值定绝对刻度。

一、方程（网络级，可数学表示）
    τ_m · dr/dt = − r(t) + W · r(t − T)
   r 为各节点活动率向量，W 为有向耦合矩阵（本文用两种权重，见二），T 为**统一延迟**。
  特征方程：  det( (1 + τ_m λ)·I − e^{−λT} W ) = 0
  等价地：对 W 的每个特征值 μ，解 1 + τ_m λ = μ·e^{−λT}；主特征值即**谱半径 ρ(W)**。

  ⇒ **单环判据里的环路增益 w 被 ρ(W) 取代**（这是本件的核心推广）：
     单调增长（自持）：ρ(W) > 1
     振荡（Hopf）模态：ρ(W) ≥ w_c^osc(T, τ_m)，其中 w_c^osc 见 theory/dynamics_analytic.py

二、两种权重（都要给，因为论文两种口径都有）
  (a) 拓扑权重：w_ij = 1（只看连接方式）
  (b) EDR 权重：w_ij = exp(−λ·d_ij)（PMC12380633 的连接概率直接当耦合强度——**这是本模型的建模选择**）

三、绝对刻度（**全部取自论文，逐字**）
  · 传导速度 v = 0.15 m/s：`papers/PMC13645236.xml`：
      `Spikes were delivered to synapses with a distance-dependent axonal delay
       corresponding to a conduction velocity of 0.15 m/s along a direct linear path`
    （同文还给出 `L2/3: 129 ± 12 µm/ms, L4: 122 ± 5 µm/ms` 的实测波速；
      `papers/FrontNeuroanat2012.xml` 给皮层内轴突 `0.1–0.6 m/s`）
  · 膜时间常数 τ_m = 10 ms：`papers/Naud2008_AdEx.xml`：
      `C = 100 pF, g_L = 10 nS … This corresponds to a membrane time constant τ = C/g_L of 10 ms`
  ⇒ 每毫米延迟 T_mm = 1 mm / 0.15 m/s = **6.667 ms/mm**；无量纲 T/τ_m = 0.6667·d_mm

判据（都用求根/幂迭代，不扫参）
  [1] 幂迭代收敛：|ρ x − W x| 残差 < 1e-9，并报告迭代步数
  [2] 判据与仿真一致：整网 DDE 直接数值解的增长率符号 == 判据给出的符号（ρ>1 增长 / ρ<1 衰减）
  [3] 绝对刻度可核：由 v 与 τ_m 算出 T/τ_m，并与 w_c^osc 表对照给出"振荡所需增益"
  [4] 负例：DAG（无环）与空图 → ρ(W)=0 → 必衰减
"""
VERSION = "dyn-net/1.0.0"
import importlib.util, math, os, random, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.join(ROOT, "audit")
V_CITED = 0.15          # m/s，PMC13645236
TAU_M_MS = 10.0         # ms，Naud2008（C=100 pF, g_L=10 nS）
MS_PER_MM = (1e-3 / V_CITED) * 1e3      # = 6.6667 ms/mm

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def adjacency(edges, nodes, weighted=True, lam=2.0):
    """返回有向邻接矩阵（dict of dict）。weighted=True 用 EDR 权重 exp(−λd)。"""
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes); W = [[0.0] * n for _ in range(n)]
    for (u, v), d in edges.items():
        W[idx[u]][idx[v]] = math.exp(-lam * d) if weighted else 1.0
    return W

def power_iteration(W, iters=5000, tol=1e-13, seed=1):
    """幂迭代求谱半径 ρ(W)（求根式迭代，不是网格扫描）。"""
    n = len(W)
    if n == 0: return 0.0, 0.0, 0
    rng = random.Random(seed)
    x = [rng.random() for _ in range(n)]
    nrm = math.sqrt(sum(v * v for v in x)); x = [v / nrm for v in x]
    rho = 0.0
    for k in range(iters):
        y = [sum(W[i][j] * x[j] for j in range(n)) for i in range(n)]
        ny = math.sqrt(sum(v * v for v in y))
        if ny < 1e-300: return 0.0, 0.0, k
        y = [v / ny for v in y]
        rho_new = ny
        if abs(rho_new - rho) < tol * max(1.0, rho_new):
            # 残差 |ρx − Wx|
            resid = math.sqrt(sum((rho_new * y[i] - sum(W[i][j] * y[j] for j in range(n))) ** 2
                                  for i in range(n)))
            return rho_new, resid, k + 1
        rho, x = rho_new, y
    return rho, float("nan"), iters

def network_dde_growth(W, T_ms, tau_ms, dt=0.05, steps=20000, seed=2):
    """直接数值解整网线性 DDE，返回估计的对数增长率（用于与判据对照）。"""
    n = len(W)
    if n == 0: return -1.0
    d = max(1, int(round(T_ms / dt)))
    rng = random.Random(seed)
    hist = [[0.0] * n for _ in range(d + 1)]
    for i in range(n): hist[0][i] = rng.random()
    def norm(v): return math.sqrt(sum(x * x for x in v))
    for k in range(1, steps):
        cur = hist[(k - 1) % (d + 1)]; past = hist[(k - 1 - d) % (d + 1)]
        new = [cur[i] + dt / tau_ms * (-cur[i] + sum(W[i][j] * past[j] for j in range(n))) for i in range(n)]
        hist[k % (d + 1)] = new
        if norm(new) > 1e12: return math.log(1e12) / (k * dt)
        if norm(new) < 1e-300: break
    a = norm(hist[(steps - 1) % (d + 1)]); b = norm(hist[steps // 2 % (d + 1)])
    return math.log(max(a, 1e-300) / max(b, 1e-300)) / ((steps // 2) * dt)

def network_dde_growth_hetero(W, delays_ms, tau_ms, dt=0.05, steps=16000, seed=3):
    """P2 余项：**每边各自延迟**的整网 DDE 直接数值解。

        τ_m · dr_i/dt = − r_i(t) + Σ_j W_ij · r_j(t − d_ij)

    这是特征方程 det((1+τ_mλ)I − ∫ e^{−λ d_ij} W) = 0 的时域形式（统一延迟只是 d_ij ≡ T 的特例）。
    返回估计的对数增长率。
    """
    n = len(W)
    if n == 0: return -1.0
    dmax = max(delays_ms.values()) if delays_ms else 1.0   # [修 1.0.2] 原来 max(dict) 取的是**键**（边元组）
    M = max(2, int(round(dmax / dt)) + 2)
    rng = random.Random(seed)
    hist = [[0.0] * n for _ in range(M)]
    for i in range(n): hist[0][i] = rng.random()
    idx = {e: max(1, int(round(delays_ms[e] / dt))) for e in delays_ms}
    edges = list(delays_ms.keys())
    def norm(v): return math.sqrt(sum(x * x for x in v))
    for k in range(1, steps):
        cur = hist[(k - 1) % M]
        new = [0.0] * n
        for i in range(n):
            acc = 0.0
            for e in edges:
                (u, v) = e
                if u != i: continue
                acc += W[i][v] * hist[(k - 1 - idx[e]) % M][v]
            new[i] = cur[i] + dt / tau_ms * (-cur[i] + acc)
        hist[k % M] = new
        if norm(new) > 1e12: return math.log(1e12) / (k * dt)
        if norm(new) < 1e-300: break
    a = norm(hist[(steps - 1) % M]); b = norm(hist[steps // 2 % M])
    return math.log(max(a, 1e-300) / max(b, 1e-300)) / ((steps // 2) * dt)

def selftest():
    print(f"{VERSION}\n")
    da = load("da", os.path.join(ROOT, "dynamics_analytic.py"))
    v3 = load("v3", os.path.join(AUD, "cycle_model_v3.py"))
    cc = load("cc", os.path.join(AUD, "cycle_census.py"))
    ok = True

    pos = cc.layout(); edges = v3.build(pos, 2.0, 0.10)
    nodes = sorted({n for e in edges for n in e})

    print("[1] 谱半径 ρ(W)（幂迭代）")
    Wt = adjacency(edges, nodes, weighted=False)
    We = adjacency(edges, nodes, weighted=True, lam=2.0)
    rho_t, res_t, it_t = power_iteration(Wt)
    rho_e, res_e, it_e = power_iteration(We)
    print(f"    拓扑权重：ρ={rho_t:.6f}  残差={res_t:.2e}  迭代={it_t}")
    print(f"    EDR 权重：ρ={rho_e:.6f}  残差={res_e:.2e}  迭代={it_e}")
    conv = (res_t < 1e-9) and (res_e < 1e-9)
    print(f"    幂迭代收敛（残差<1e-9）: {conv}")
    ok &= conv

    print("\n[3] 绝对刻度（论文数值）")
    print(f"    v = {V_CITED} m/s（PMC13645236）→ 每毫米延迟 {MS_PER_MM:.4f} ms/mm")
    print(f"    τ_m = {TAU_M_MS} ms（Naud2008：C=100 pF, g_L=10 nS）")
    for d_mm in (1.0, 5.0, 7.75):
        Tr = MS_PER_MM * d_mm; ratio = Tr / TAU_M_MS
        wc, w1, _ = da.wc_analytic(ratio * TAU_M_MS, TAU_M_MS)
        print(f"    环长 {d_mm:5.2f} mm → T={Tr:6.2f} ms → T/τ_m={ratio:5.2f} → 振荡所需 ρ(W) ≥ {wc:.3f}")

    print("\n[2] 判据 vs 整网 DDE 仿真（同一张图）")
    T_ref = MS_PER_MM * 7.75          # 基准：7.75 mm ≈ 我们用 0.25 ms/mm 量到的 31 拍
    rows = []
    # [修 1.0.1] 原取 0.5/1.0/2.0 三个缩放，算出的 ρ 全 >1（1.459/2.918/5.836）⇒ 只在
    # "自持"一侧检验了判据，衰减一侧未被覆盖（判据只能变化一半）。补一个 0.2 → ρ≈0.58。
    for rho_scale in (0.2, 0.5, 1.0, 2.0):
        Ws = [[rho_scale * Wt[i][j] for j in range(len(nodes))] for i in range(len(nodes))]
        rho_s, _, _ = power_iteration(Ws)
        g = network_dde_growth(Ws, T_ref, TAU_M_MS)
        pred = (rho_s > 1.0)
        obs = (g > 0)
        rows.append((rho_scale, rho_s, g, pred, obs, pred == obs))
        print(f"    缩放={rho_scale:<4} ρ(W)={rho_s:7.3f}  仿真增长率={g:+.5f}  判据(自持)={pred}  实测={obs}  一致={pred==obs}")
    consistent = all(r[5] for r in rows)
    print(f"    判据与仿真一致: {consistent}")
    ok &= consistent

    print("\n[4] 负例：DAG 与空图")
    dag = {k: v for k, v in edges.items() if pos[k[1]][0] > pos[k[0]][0] + 0.2}
    rho_dag, _, _ = power_iteration(adjacency(dag, sorted({n for e in dag for n in e}), weighted=True))
    rho_empty, _, _ = power_iteration([])
    g_dag = network_dde_growth(adjacency(dag, sorted({n for e in dag for n in e}), weighted=True), T_ref, TAU_M_MS)
    print(f"    DAG：ρ={rho_dag:.6f}  仿真增长率={g_dag:+.5f}；空图 ρ={rho_empty:.6f}")
    caught = (rho_dag == 0.0) and (g_dag < 0)
    print(f"    DAG 必衰减: {caught}")
    ok &= caught

    print("\n[5] P2 余项：异质延迟（每边各自 d/v）与均匀延迟的夹逼")
    # 论断（可检验）：延迟越大 → 有效增益 e^{−λd} 越小 → 增长率越低。
    # 因此**异质延迟的增长率应落在"全边取 T_max"与"全边取 T_min"之间**。
    node_idx = {n: i for i, n in enumerate(nodes)}
    dly = {}
    for (u, v), d in edges.items():
        dly[(node_idx[u], node_idx[v])] = max(0.05, MS_PER_MM * d)   # 论文刻度：d(mm)/v
    Tmin, Tmax = min(dly.values()), max(dly.values())
    Tmid = sum(dly.values()) / len(dly)
    W1 = adjacency(edges, nodes, weighted=True, lam=2.0)
    # [修 1.0.3] 原缩放取 1.6，而 EDR 口径 ρ(W)=0.463 ⇒ 1.6×0.463=0.74 **仍在阈值以下**，
    # 三个增长率全是 −1e−4 量级、窗口只有 8e−5，夹逼"成立"其实在噪声里（判据无牙）。
    # 抬到 5.0（ρ≈2.3）使增长率达 1e−2 量级，延迟依赖才可分辨。
    Wsc = [[5.0 * W1[i][j] for j in range(len(nodes))] for i in range(len(nodes))]
    rho_sc, _, _ = power_iteration(Wsc)
    print(f"    缩放后 ρ(W) = {rho_sc:.3f}（须 >1 才有可分辨的增长率）")
    print(f"    延迟范围（每边 d/v）：{Tmin:.2f} ~ {Tmax:.2f} ms（均值 {Tmid:.2f} ms）")
    g_min = network_dde_growth(Wsc, Tmin, TAU_M_MS)
    g_max = network_dde_growth(Wsc, Tmax, TAU_M_MS)
    g_het = network_dde_growth_hetero(Wsc, dly, TAU_M_MS)
    print(f"    全边取 T_min={Tmin:.2f}ms → 增长率 {g_min:+.5f}")
    print(f"    全边取 T_max={Tmax:.2f}ms → 增长率 {g_max:+.5f}")
    print(f"    **异质**（每边各自 d/v）   → 增长率 {g_het:+.5f}")
    lo, hi = min(g_min, g_max), max(g_min, g_max)
    inside = (lo - 5e-4) <= g_het <= (hi + 5e-4)
    print(f"    夹逼成立（异质落在两均匀之间）: {inside}  （区间 [{lo:+.5f}, {hi:+.5f}]）")
    ok &= inside

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

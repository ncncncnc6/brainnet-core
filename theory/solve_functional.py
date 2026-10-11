#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/solve_functional.py —— 解泛函：从平衡因子到指数距离律的 λ（theory-solve-functional/1.0.0）

本件**独立于源码树**（放在 theory/ 下），不改动 brainnet-core 的任何实现件。

一、泛函的逐字出处（papers/ 目录内，可核对）
  · 形状（平衡因子）：PMC11790248（eLife 2025，树突生长）
      “The balance in the model is modulated by a balancing factor ( b f ). If b f is zero,
       dendritic trees minimize the cable only, and if b f is one, equal weight is given to
       minimizing conduction delays and total dendritic length.”
    ⇒ C = Σ_edges 长度 + bf · Σ_edges 路径长度（bf ∈ [0,1]；bf=0 纯电缆，bf=1 电缆与延迟等权）
  · 权衡的两侧与实测效率：PMC2837396（PLoS Comput Biol 2010，轴突乔木在材料与传导延迟之间权衡）
      给出最小长度图用 Euclidean Steiner minimal tree 度量，并报告 spiny 细胞的效率
      ε_spiny = 0.86 ± 0.04（即真实乔木**并非**最优，而是接近最优）。
  · 距离衰减的形式：PMC12380633 式(1) ℙ(J=1)=exp(−λ_decay·d)，其 λ 的验证区间 1.2–3.2 mm⁻¹。

二、解泛函（本件自己做的推导，逐步可核）
  把"选哪些连接"松弛成距离上的分布 p(d)，则代价约束下的最优分布由**自由能**给出：
      F[p] = Σ_i p_i · c_i  −  T · H(p),      c_i = (1 + bf) · d_i,   H(p) = −Σ p_i ln p_i
  其中 c_i 是单位连接在距离 d_i 处的代价：长度项 d_i 加上经 bf 加权的路径（延迟）项 bf·d_i；
  T 是"无序度尺度"（与距离同量纲，mm）。对 Σp=1 求极小：
      指数族 ⇒  p(d) ∝ exp(−(1+bf)·d / T)
  即**指数距离律，且 λ = (1 + bf) / T**。这就是泛函的 Lagrange 乘子与 λ 的对应关系。

三、判据
  [1] 数值解：用**镜像下降**（不是梯度下降；见 solve 的注释）在单纯形上最小化 F[p]，
      拟合 log p vs d 的斜率 → 应等于 (1+bf)/T
  [2] bf 扫描：T 固定时 λ 随 bf 线性上升（λ = (1+bf)/T）
  [3] 反演：由文献 λ=2.0 mm⁻¹ 反解 (bf, T) 的可行组合
  [4] 负例：T=1e-6（**近似**去掉熵项，不是 T 严格取 0）→ 分布退化到最短距离处，
      斜率拟合失效——判据必须察觉。
  [5] **独立的变分核对（P5）**：让 `free_energy` 真正参与判据 ——
      在解出的 p̂ 上比较 F[p̂] 与 F[解析解]，并做**扰动最优性**检验：
      任何可行扰动都必须使 F 上升；再以"故意用错 λ"作负例。这是与求解器更新式**无关**的检验，
      因此不再是构造性自洽（D 组指出：镜像下降的不动点本就与解析解同构，0.00% 因此不算证据）。
"""
VERSION = "theory-solve-functional/1.0.0"
import math, sys

D_LO, D_HI, NB = 0.01, 2.20, 600

def grid():
    h = (D_HI - D_LO) / NB
    return [D_LO + (i + 0.5) * h for i in range(NB)], h

def free_energy(p, ds, bf, T, h=None):
    """F[p] = Σ p_i·c_i − T·H(p)，c_i=(1+bf)·d_i，H(p) = −Σ p_i ln p_i。

    [修 1.0.1 / P5] 原实现把**两项都乘 h**（−T·H 写成 +T·Σp ln p·h），量纲与符号都不自洽，
    且从未被任何判据调用（D 组 41-46：0.00% 因此是构造性自洽而非独立验证）。
    现按定义写：代价项一阶、熵项带负号，h 参数保留但**不再使用**（仅为签名稳定）。
    """
    F = sum(pi * (1 + bf) * d for pi, d in zip(p, ds))
    H = -sum(pi * math.log(pi) for pi in p if pi > 0)
    return F - T * H

def solve(bf, T, iters=3000, ds=None, h=None):
    """在单纯形上用**镜像下降**最小化 F[p]（自写，无外部依赖）。

    [1.0.1] 首版用 z 空间梯度下降并把步长乘了 bin 宽 h(≈0.0037) → 实际步长 ~1e-3，
    3000 轮后分布仍近似均匀（max p≈0.002），拟合出的 λ≈0.01，与解析值差 99%。
    改为镜像下降：p ← p·exp(−τ·g) / Z，τ = 1/(1+T)；固定点即指数族 p∝exp(−c/T)。
    """
    if ds is None: ds, h = grid()
    n = len(ds)
    p = [1.0 / n] * n
    tau = 1.0 / (1.0 + T)
    for _ in range(iters):
        g = [(1 + bf) * d + T * (math.log(pi + 1e-300) + 1.0) for d, pi in zip(ds, p)]
        w = [pi * math.exp(-tau * gi) for pi, gi in zip(p, g)]
        S = sum(w) or 1.0
        p = [wi / S for wi in w]
    return p

def fit_lambda(p, ds):
    """对 log p 与 d 做最小二乘（只用 p 明显的区段，避免数值下溢端）。
    [1.0.2] T→0 时分布塌到单个 bin，可用点只有 1 个 → 原来的 num/den 直接除零。
    改为返回 None 表示"不可拟合"（塌缩本身就是要检测的现象）。"""
    xs, ys = [], []
    pmax = max(p)
    for d, pi in zip(ds, p):
        if pi > 1e-9 and pi > pmax * 1e-4:
            xs.append(d); ys.append(math.log(pi))
    if len(xs) < 2: return None
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den <= 0: return None
    return -(num / den)                                    # 斜率取负 = λ

def selftest():
    print(f"{VERSION}\n")
    ds, h = grid(); ok = True

    print("[1] 数值解 vs 解析解 λ = (1+bf)/T")
    for bf, T in ((0.0, 1.0), (0.5, 1.0), (1.0, 1.0), (1.0, 0.5)):
        p = solve(bf, T, ds=ds, h=h)
        lam_num = fit_lambda(p, ds); lam_ana = (1 + bf) / T
        err = abs(lam_num - lam_ana) / lam_ana
        good = err < 0.05; ok &= good
        print(f"    bf={bf:<4} T={T:<4} λ数值={lam_num:6.3f}  λ解析={lam_ana:6.3f}  相对误差={err*100:5.2f}%  {'OK' if good else 'FAIL'}")

    print("\n[2] bf 扫描（T=1.0 mm 固定）：λ 随 bf 线性上升")
    lams = []
    for bf in (0.0, 0.25, 0.5, 0.75, 1.0):
        lam = fit_lambda(solve(bf, 1.0, ds=ds, h=h), ds); lams.append(lam)
        print(f"    bf={bf:<5} λ={lam:6.3f}   解析={(1+bf)/1.0:6.3f}")
    mono = all(lams[i] < lams[i + 1] for i in range(len(lams) - 1)); ok &= mono
    print(f"    单调递增: {mono}")

    print("\n[3] 反演：文献 λ=2.0 mm⁻¹ 对应的 (bf, T) —— 真判据（含失败分支）")
    target = 2.0; inv_ok = True
    for bf in (0.0, 0.5, 1.0):
        T = (1 + bf) / target
        lam_fit = fit_lambda(solve(bf, T, ds=ds, h=h), ds)
        err = abs(lam_fit - target) / target
        good = err < 0.05
        inv_ok &= good
        print(f"    bf={bf:<4} ⇒ T=(1+bf)/λ={T:.3f} mm → 解出 λ={lam_fit:.3f} vs 目标 {target}  相对误差={err*100:5.2f}%  {'OK' if good else 'FAIL'}")
    # 负例：故意用错 bf，反演必须失败
    lam_wrong = fit_lambda(solve(1.5, (1 + 0.0) / target, ds=ds, h=h), ds)
    teeth3 = abs(lam_wrong - target) / target > 0.05
    inv_ok &= teeth3
    print(f"    负例：bf 故意用错（1.5 配 bf=0 的 T）→ λ={lam_wrong:.3f}，偏离超 5%: {teeth3}")
    ok &= inv_ok

    print("\n[4] 负例：T=0（删掉熵项）→ 分布必须塌缩，斜率不可拟合")
    p0 = solve(1.0, 1e-6, ds=ds, h=h)
    lam0 = fit_lambda(p0, ds)
    collapsed = (max(p0) > 0.5) or (lam0 is None)
    caught = collapsed; ok &= caught
    ltxt = "不可拟合（单点）" if lam0 is None else f"{lam0:.1f}"
    print(f"    T≈0：λ数值={ltxt}  最大概率质量={max(p0):.4f} → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n[5] 独立变分核对（P5）：free_energy 首次参与判据")
    import random as _r
    bf0, T0 = 1.0, 1.0
    lam0_ = (1 + bf0) / T0
    p_hat = solve(bf0, T0, ds=ds, h=h)
    w = [math.exp(-lam0_ * d) for d in ds]; Z = sum(w)
    p_ana = [wi / Z for wi in w]
    F_hat = free_energy(p_hat, ds, bf0, T0)
    F_ana = free_energy(p_ana, ds, bf0, T0)
    rel = abs(F_hat - F_ana) / abs(F_ana)
    print(f"    F[p̂]={F_hat:.6f}  F[p_解析]={F_ana:.6f}  相对差={rel*100:.4f}%")
    # 扰动最优性：任何可行扰动都必须使 F 上升
    rng = _r.Random(3); worse = 0; trials = 40
    for _ in range(trials):
        q = list(p_hat)
        for _k in range(5):
            i, j = rng.randrange(len(q)), rng.randrange(len(q))
            eps = 1e-4 * rng.random() * q[i]
            q[i] -= eps; q[j] += eps
        S = sum(q); q = [x / S for x in q]
        if free_energy(q, ds, bf0, T0) >= F_hat: worse += 1
    pert_ok = (worse == trials)
    print(f"    扰动最优性：{worse}/{trials} 个随机可行扰动的 F 都不降低 → {pert_ok}")
    # 负例：故意用错 λ 解析解，F 必须更高
    lam_bad = 2 * lam0_
    wb = [math.exp(-lam_bad * d) for d in ds]; Zb = sum(wb)
    F_bad = free_energy([x / Zb for x in wb], ds, bf0, T0)
    teeth5 = (F_bad > F_hat)
    print(f"    负例：λ 取 2 倍（{lam_bad:.2f}）→ F={F_bad:.6f} > F[p̂]={F_hat:.6f}: {teeth5}")
    variational_ok = (rel < 1e-3) and pert_ok and teeth5
    print(f"    → 独立变分核对: {variational_ok}（该检验**不使用**求解器的更新式）")
    ok &= variational_ok

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

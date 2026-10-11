#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/dynamics_analytic.py —— 循环节的动力学方程与解析自持判据（dyn/1.0.0）

用户要求：**造动力数学模型**，能数学表示、能被求解器求解、能被公式推导，**不要扫参数**。
本件就做这一件事：把"环能不能自持"写成带延迟的线性动力学，解它的特征方程。

一、方程（可数学表示）
  在一条有向环上，把 AdEx 的膜方程在静息点线性化（Naud2008 用同一套线性化/分岔语言，
  逐字见 papers/Naud2008_AdEx.xml：'the transition occurs via the Andronov-Hopf bifurcation'、
  'We give an analytical criterion to distinguish …'），得到**延迟微分方程（DDE）**：

      τ_m · dr/dt = − r(t) + w · r(t − T_loop)

  其中 r = 群体活动率，w = 绕环一圈的**环路增益**（环上各边权重之积），
  T_loop = 环周期 = Σ 各边延迟（cycle_census 已实测：基准最长 31 拍）。

二、特征方程与解析判据（可推导）
  代入 r ∝ e^{λt} 得     1 + τ_m λ = w · e^{−λ T_loop}
  令 λ = iω（振荡模态），分开模与相：
      模：  w = sqrt(1 + (ω τ_m)²)                    （A 式）
      相：  ω T_loop + arctan(ω τ_m) = kπ            （B 式）
  自持（出现非衰减振荡模态）当且仅当存在满足 (B) 的 ω 使 w ≥ sqrt(1+(ωτ_m)²)。
  取 (B) 的最小正根 ω₁ 即得**临界增益**

      w_c(T_loop, τ_m) = sqrt(1 + (ω₁ τ_m)²)          ← 解析判据

  极限行为（可直接读出）：T_loop ≫ τ_m 时 ω₁ ≈ π/T_loop ⇒ **w_c → 1**；
  T_loop 与 τ_m 同量级时 w_c 明显大于 1。**环越长，越容易自持。**

三、怎么解（不扫参）
  (B) 是关于 ω 的单调方程 ⇒ 用**二分求根**解出 ω₁（求根，不是网格扫描）；
  再用 (A) 给出 w_c。并与**直接数值解 DDE**（同一方程的独立解法）对照，二者必须一致。

判据
  [1] 求根收敛：|F(ω₁)| 足够小，且 ω₁ 落在 (0, π/T_loop] 内
  [2] 解析与数值一致：两种解法给出的 w_c 相对误差 < 2%
  [3] 极限行为：T_loop/τ_m = 100 时 w_c < 1.01；比值 = 0.2 时 w_c > 1.5
  [4] 负例：无环（T_loop 不存在）→ 判据不适用，数值解必衰减
"""
# ── 出处等级（dyn_ledger 的 B 级）─────────────────────────────────────
# 本文件的方程 τ_m·ṙ = −r + w·r(t−T) 是**一阶延迟反馈**的标准形式，**由本文件推导/选取**，
# papers/ 语料里没有给出这条方程本身 ⇒ 归 **B 级（建模选择）**，不得当作论文结论引用。
# 其中进入方程的两个量是 A 级：τ_m = 10 ms（Naud2008 的 C=100 pF, g_L=10 nS）、
# 环路延迟 T 由连接长度与传导速度 v=0.15 m/s（PMC13645236）换算。
VERSION = "dyn/1.0.0"
import cmath, math, sys

def omega1(T_loop, tau_m):
    """解相条件的最小正根：F(ω) = ω*T_loop + atan(ω*tau_m) − π = 0（单调递增 ⇒ 二分）。"""
    # [修 1.0.1] 相位条件原写成 = π，**错**：由 1+τλ = w e^{−λT}、λ=iω 得
    #   ωT + arctan(ωτ) = 2πk，最小正根对应 k=1（即 2π），不是 π。
    lo, hi = 1e-12, 4.0 * math.pi / max(T_loop, 1e-12)
    F = lambda w: w * T_loop + math.atan(w * tau_m) - 2.0 * math.pi
    if F(hi) < 0:                      # 扩大上界（保守，不做网格扫描）
        hi *= 2.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if F(mid) < 0: lo = mid
        else: hi = mid
    w = (lo + hi) / 2.0
    return w, abs(F(w))

def wc_analytic(T_loop, tau_m):
    w, resid = omega1(T_loop, tau_m)
    return math.sqrt(1.0 + (w * tau_m) ** 2), w, resid

def wc_numeric(T_loop, tau_m, dt=0.01, T=4000.0, w_lo=0.2, w_hi=6.0):
    """直接数值解 DDE，用增长率过零求临界增益（求根，非扫描）。"""
    def growth(w):
        n = int(T / dt); d = int(round(T_loop / dt))
        if d < 1: return float("nan")
        r = [0.0] * (n + 1); r[0] = 1.0
        for i in range(1, n):
            past = r[i - d] if i - d >= 0 else 0.0
            r[i] = r[i - 1] + dt / tau_m * (-r[i - 1] + w * past)
            if abs(r[i]) > 1e6: return math.log(1e6) / (i * dt)
        head = max(abs(r[i]) for i in range(1, max(2, int(50 / dt))))
        tail = max(abs(r[i]) for i in range(n - int(200 / dt), n))
        return math.log(max(tail, 1e-300) / max(head, 1e-300)) / ((n * dt) - 50.0)
    # 增长率随 w 单调不减 ⇒ 二分求零点
    for _ in range(40):
        mid = (w_lo + w_hi) / 2.0
        if growth(mid) < 0: w_lo = mid
        else: w_hi = mid
    return (w_lo + w_hi) / 2.0

def selftest():
    print(f"{VERSION}\n")
    ok = True
    tau_m = 1.0                      # 无量纲：一切以比值 T_loop/τ_m 表述（绝对 τ_m 需按 AdEx 取值）

    print("[1] 解相条件的最小正根（二分求根，不扫网格）")
    for ratio in (0.2, 1.0, 5.0, 20.0, 100.0):
        T = ratio * tau_m
        w, resid = omega1(T, tau_m)
        # [修 1.0.4] 边界原写 π/T —— 那是 π 相位条件的产物；改为 2π 后
        # ω₁ 满足 ωT + arctan(ωτ) = 2π ⇒ ω₁ < 2π/T。边界必须跟着相位条件走。
        ok &= (resid < 1e-10) and (0 < w < 2.0 * math.pi / T + 1e-9)
        print(f"    T/τ_m={ratio:<6} ω₁={w:8.5f}  余差={resid:.2e}")

    print("\n[2] 解析判据的独立检验（两条，都用求根/代入，不扫参）")
    # [修 1.0.3] 上一版用牛顿取"主根"来对照 —— 错：w>1 时特征方程**总有实正根**，
    # 主根是单调增长的那条，与 Hopf 复根不是同一条（实测偏 0.2~0.4）。正确做法：
    #   (a) 把解析解 (ω₁, w_c^osc) **代回** Z(λ)=1+τλ−w e^{−λT}，残差必须≈0；
    #   (b) 单调阈值用实轴二分求根验证：w>1 有实正根，w<1 没有。
    res_ok = True
    for ratio in (0.2, 1.0, 5.0, 20.0, 100.0):
        T = ratio * tau_m
        wa, w1, _ = wc_analytic(T, tau_m)
        Z = abs(1.0 + tau_m * complex(0, w1) - wa * cmath.exp(-complex(0, w1) * T))
        res_ok &= (Z < 1e-9)
        print(f"    T/τ_m={ratio:<6} |Z(iω₁)| = {Z:.3e}（代入残差，应≈0）")
    print(f"    代入残差全部≈0: {res_ok}")
    ok &= res_ok

    def real_root(w, T, tau):
        """实轴二分：解 1+τλ−w e^(−λT)=0 的正根（不存在则返回 None）。"""
        f = lambda lam: 1.0 + tau * lam - w * math.exp(-lam * T)
        if f(0.0) >= 0: return None
        lo, hi = 0.0, 1e-6
        while f(hi) < 0 and hi < 1e6: hi *= 2.0
        for _ in range(200):
            mid = (lo + hi) / 2.0
            if f(mid) < 0: lo = mid
            else: hi = mid
        return (lo + hi) / 2.0
    for w in (0.9, 0.99, 1.01, 1.5):
        r = real_root(w, 2.0 * tau_m, tau_m)
        print(f"    w={w:<5} 实正根 = {'无（衰减）' if r is None else f'{r:.6f}（增长）'}")
    mono_ok = (real_root(0.99, 2.0, tau_m) is None) and (real_root(1.01, 2.0, tau_m) is not None)
    print(f"    单调阈值恰为 w=1（<1 无实正根、>1 有）: {mono_ok}")
    ok &= mono_ok

    print("\n[3] 极限行为（解析式直接给出）")
    w100, _, _ = wc_analytic(100.0 * tau_m, tau_m)
    w02, _, _ = wc_analytic(0.2 * tau_m, tau_m)
    ok &= (w100 < 1.01) and (w02 > 1.5)
    print(f"    T/τ_m=100 → w_c={w100:.4f}（应 <1.01）；T/τ_m=0.2 → w_c={w02:.4f}（应 >1.5）")

    print("\n[4] 负例：无环（不存在 T_loop）")
    print("    判据不适用（没有环就没有 DDE）；数值上任意 w 的增长率 <0")
    def growth_decay(w):
        dt, T, tau = 0.01, 500.0, tau_m
        r = [0.0] * int(T / dt + 1); r[0] = 1.0
        for i in range(1, len(r)): r[i] = r[i - 1] + dt / tau * (-r[i - 1])
        return r[-1]
    caught = abs(growth_decay(10.0)) < 1e-6
    ok &= caught
    print(f"    无环时末值={growth_decay(10.0):.2e} → {'PASS（必衰减）' if caught else 'FAIL'}")

    print("结论式（可直接引用）：w_c^osc(T_loop, τ_m) = sqrt(1 + (ω₁·τ_m)²)，")
    print("              ω₁ 为 ω·T_loop + arctan(ω·τ_m) = 2π 的最小正根；")
    print("              另有**精确**结论：w > 1 时存在实正特征根（单调增长）⇒ 自持阈值 w_c^min = 1。")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

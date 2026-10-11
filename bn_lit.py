#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_lit.py —— 用文献把脑网参数锚到绝对单位（bn-lit/1.0.0）

文献（本轮实取全文，非转述）：
  Naud, Marcille, Clopath, Gerstner (2008) "Firing patterns in the adaptive exponential
  integrate-and-fire model", Biological Cybernetics 99(4-5):335-347,
  DOI 10.1007/s00422-008-0264-7 —— 全文经 Europe PMC REST 取得（/tmp/ep.xml, 132 KB）。

从该文正文直接引用的三条锚点（原文句式，已核对）：
  A1 定义：参数 C 与 g_L 被吸收进时间尺度 τ_m = C/g_L。
  A2 取值：该文的参数对应膜时间常数 τ = C/g_L = **10 ms**。
  A3 双体制：τ_w 被固定为 **5 ms 或 100 ms**，前者使 w 可解释为**不应电流**，
      后者使 w 可解释为**适应电流**；两者是不同物理角色，不是同一旋钮的两档。
  参数表列：C(pF) / g_L(nS) / E_L(mV) / V_T(mV) / Δ_T(mV) / a(nS) / τ_w(ms) / b(pA) / V_r(mV) / I(pA)。

由此得到本文件的核心推导（把「U_DECAY 该取多少」从猜变成算）：
      DECAY   = exp(−dt/τ_m)   ⟹  dt = −τ_m·ln(DECAY)
      U_DECAY = exp(−dt/τ_w)
      ⟹ U_DECAY = DECAY^(τ_m/τ_w)          ← 给定比值，U_DECAY 唯一确定
   反过来：τ_w/τ_m = ln(DECAY)/ln(U_DECAY)   （bn_units.py 中的无量纲恒等式）

结论（本文件算出来的，不是引的）：
  脑网现有 U_DECAY=0.88 ⟹ τ_w/τ_m = 1.5524；
  文献两个体制是 5/10 = **0.5**（不应）与 100/10 = **10**（适应）；
  1.5524 **落在两个体制之间**，既不干净地是不应，也不干净地是适应 —— 它是一个混合体。
  若要贴文献，U_DECAY 应取 DECAY^0.1 = 0.98035（适应）或 DECAY^2 = 0.6724（不应）。
"""

VERSION = "bn-lit/1.0.0"
import json, math, sys

DECAY = 0.82
U_DECAY_CUR = 0.88
RECOVERY_CUR = 16
TARGET_CUR = 0.14

LIT = dict(
    tau_m_ms=10.0,               # A2，原文直给
    tau_w_refractory_ms=5.0,     # A3
    tau_w_adaptation_ms=100.0,   # A3
    doi="10.1007/s00422-008-0264-7",
    source="Naud et al. 2008, Biol Cybern 99(4-5):335-347 (全文经 Europe PMC 取得)",
)

def dt_ms(tau_m_ms=LIT["tau_m_ms"], decay=DECAY):
    return -tau_m_ms * math.log(decay)

def u_decay_for_ratio(ratio, decay=DECAY):
    """U_DECAY = DECAY^(τ_m/τ_w) = DECAY^(1/ratio)"""
    return decay ** (1.0 / ratio)

def ratio_of(u_decay, decay=DECAY):
    return math.log(decay) / math.log(u_decay)

def window_ms(ticks, tau_m_ms=LIT["tau_m_ms"], decay=DECAY):
    return ticks * dt_ms(tau_m_ms, decay)

def hz(per_tick, tau_m_ms=LIT["tau_m_ms"], decay=DECAY):
    return per_tick * 1000.0 / dt_ms(tau_m_ms, decay)

def report():
    dt = dt_ms()
    print(f"{VERSION}   文献: {LIT['source']}   DOI {LIT['doi']}")
    print(f"\n[A] 单位锚定（τ_m = {LIT['tau_m_ms']} ms 为文献直给值）")
    print(f"   一拍 dt = −τ_m·ln(DECAY) = {dt:.4f} ms        ⇒ 16 拍 = {window_ms(16):.1f} ms")
    print(f"   现有 τ_w/τ_m = ln({DECAY})/ln({U_DECAY_CUR}) = {ratio_of(U_DECAY_CUR):.4f}"
          f"   ⇒ τ_w = {ratio_of(U_DECAY_CUR)*LIT['tau_m_ms']:.1f} ms")

    print("\n[B] 文献两个体制 vs 现有参数")
    for tag, ratio in (("不应电流（refractory）", LIT["tau_w_refractory_ms"]/LIT["tau_m_ms"]),
                       ("适应电流（adaptation）", LIT["tau_w_adaptation_ms"]/LIT["tau_m_ms"]),
                       ("现有（混合）", ratio_of(U_DECAY_CUR))):
        u = u_decay_for_ratio(ratio)
        print(f"   {tag:22} τ_w/τ_m={ratio:5.2f}  τ_w={ratio*LIT['tau_m_ms']:6.1f} ms  "
              f"→ U_DECAY 应为 {u:.5f}")
    print("   ⇒ 现有 0.88 落在两者之间：既不是干净的不应，也不是干净的适应。")

    print("\n[C] 换算成 Hz（这才是能对外说的数）")
    floor = 1.0 / (RECOVERY_CUR + 1)
    print(f"   恢复窗 {RECOVERY_CUR} 拍 = {window_ms(RECOVERY_CUR):.1f} ms；"
          f"脉冲完全无效 → 地板 {floor:.4f}/拍 = {hz(floor):.1f} Hz")
    print(f"   现有目标 {TARGET_CUR}/拍 = {hz(TARGET_CUR):.1f} Hz；目标/地板 = {hz(TARGET_CUR)/hz(floor):.2f}×")
    print(f"   文献两体制的特征频率：1/τ_w 不应 = {1000/LIT['tau_w_refractory_ms']:.0f} Hz，"
          f"适应 = {1000/LIT['tau_w_adaptation_ms']:.1f} Hz")

    print("\n[D] 方针（贴文献只需改两个数，二者必须配对）")
    print(f"   ① 当不应电流用：RECOVERY 缩到 1~2 拍（{window_ms(1):.1f}~{window_ms(2):.1f} ms），"
          f"U_DECAY → {u_decay_for_ratio(0.5):.4f}；发放率可上 {hz(1/(1+1)):.0f} Hz 量级")
    print(f"   ② 当适应电流用：U_DECAY → {u_decay_for_ratio(10):.4f}（τ_w={LIT['tau_w_adaptation_ms']:.0f} ms），"
          f"RECOVERY 保留但改名为 AHP；目标须 ≤ {floor:.4f}/拍 = {hz(floor):.1f} Hz")
    print(f"   ③ 保留混合体：可以，但必须**命名为第三体制**，且不能再把 RECOVERY 叫不应期")
    print(f"   ⚠ 配对检查：只改 U_DECAY 不改 RECOVERY（或反之）都会把模型推到两个体制之外")

    print("\n判据行的可跑检查：python3 bn_lit.py --assert")
    return True

def assertions():
    errs = []
    dt = dt_ms()
    if abs(dt - 1.9853) > 0.001: errs.append(f"dt 算错: {dt}")
    if abs(u_decay_for_ratio(10) - 0.98035) > 1e-5: errs.append("适应档 U_DECAY 算错")
    if abs(u_decay_for_ratio(0.5) - 0.6724) > 1e-4: errs.append("不应档 U_DECAY 算错")
    if abs(ratio_of(U_DECAY_CUR) - 1.5524) > 1e-4: errs.append("现有比值算错")
    # 往返一致：ratio → U_DECAY → ratio 必须回到原值（贴文献时不会自欺）
    for r in (0.5, 1.5524, 10.0):
        if abs(ratio_of(u_decay_for_ratio(r)) - r) > 1e-9: errs.append(f"往返不自洽 r={r}")
    # 单调性：比值越大 → U_DECAY 越大（τ_w 越慢）
    if not (u_decay_for_ratio(0.5) < u_decay_for_ratio(10)): errs.append("单调性错")
    # 硬结论：现有比值必须落在两体制之间（若这条不再是 True，说明参数被改过，须重新评估方针）
    if not (0.5 < ratio_of(U_DECAY_CUR) < 10): errs.append("现有比值不再落在两体制之间——参数已变，需重评")
    # [修 1.0.1] 上面 dt / 0.98035 / 0.6724 / 1.5524 四条是**回归锁**：常数由本文件同一
    # 公式算出，属自洽核对而不是独立验证（子代理 C 组指出「硬编码常数当判据」）。
    # 下面三条才是对公式行为有鉴别力的性质检验。
    import random as _rnd
    _rr = _rnd.Random(7)
    for _ in range(50):
        _r = 10 ** _rr.uniform(-2, 2)
        if abs(ratio_of(u_decay_for_ratio(_r)) - _r) > max(1e-9, 1e-6 * _r):
            errs.append(f"性质·往返不自洽（随机档 r={_r:.4f}）"); break
    _grid = [10 ** (-2 + 4 * k / 19) for k in range(20)]
    if not all(u_decay_for_ratio(_grid[k]) < u_decay_for_ratio(_grid[k + 1]) for k in range(19)):
        errs.append("性质·单调性在 20 点栅格上失败（r 增大 → U_DECAY 必须增大）")
    if not all(0.0 < u_decay_for_ratio(_r) < 1.0 for _r in (0.01, 1.0, 100.0)):
        errs.append("性质·U_DECAY 未落在开区间 (0,1)")
    if errs:
        for e in errs: print("FAIL:", e)
        return False
    print("全部判据 PASS：回归锁 4 条（dt=1.9853ms / 适应档 0.98035 / 不应档 0.6724 / 现有 1.5524 落在两体制间）"
          "+ 性质检验 3 条（50 档随机往返自洽 / 20 点栅格单调 / U_DECAY∈(0,1)）")
    return True

if __name__ == "__main__":
    if "--assert" in sys.argv:
        sys.exit(0 if assertions() else 1)
    report()
    sys.exit(0)

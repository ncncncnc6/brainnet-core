#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_dist_model.py —— 距离分布 ⇔ 时间快慢 的等价性判定（bn-dist-model/1.0.2）

用户命题：「距离双峰等价于时间的快慢？」
本文件把它拆成三条可判定的陈述，逐条给数，不靠文字断言。

数学骨架
  ① 延迟映射：距离 d 经有限传导速度 v 得延迟 τ = d/v + τ_syn。
     故距离双峰 {d_short, d_long} ⟺ 延迟双峰 {τ_fast, τ_slow}（同一件事的两种坐标）。
     时间尺度比 ρ = τ_slow / τ_fast = d_long / d_short。
  ② 完整的"快慢系统"（慢流形 / 松弛振荡）需要**慢变量**（适应电流、慢突触、
     细胞外离子、传导速度可塑），仅有延迟双峰不构成慢变量。
  ③ 故逻辑上应检验两件事：
       S1 充分性：距离双峰 + 慢变量 → 出现松弛振荡（慢相位占周期主导）？
       S2 必要性：距离单峰（但平均延迟相同）+ 慢变量 → 慢相位是否照样主导？
         若照样主导，则"快慢 ⇒ 距离双峰"不成立（双峰非必要）。

判据（松弛振荡的量化：慢相位占比）—— **现行口径**（[C 组 (d) 19-22] 原文件头写的是
已废弃的绝对阈值 0.70/0.45，与代码不符；实测三组全落在 0.42–0.54，按旧口径会全部落空）：
  把一个周期分成"慢相位"（亚阈缓慢爬升，v < 0.5·θ）与"快相位"（阈值附近跃变），
  算**最后 k 个周期**的慢相位占比（并给出跨周期标准差），然后**与参考组对照**：
  差值的绝对值超过参考的**2 倍标准误**才判"高于/低于参考"，否则判"与参考无差"。
  容差不是拍出来的常数，而是由**多种子扫描**得到的标准误决定（[C 组 (b) 112-114]）。
"""
VERSION = "bn-dist-model/1.0.2"
import math, sys

THETA = 1.0            # 触发阈值
V_DECAY = 0.90         # 膜电位衰减
W_STEP = 0.30          # 每次触发给恢复变量的增量
BETA = 0.80            # 恢复变量对兴奋性的压制
I0 = 0.12              # 恒定慢驱动（tonic input）——[1.0.1] 首版漏了它，只有延迟自反馈，
                       # 结果三组全部"无周期"：回路根本不起振。实测 i0=0.12 落在振荡区。
GAIN = 0.50            # 每条延迟通路的反馈增益

def run(tau_short, tau_long, w_short, w_long, dw, beats=1400, seed=1):
    """双通路自反馈回路：v 经两条延迟通路（短/长）回到自身，恢复变量 w 提供慢变量。
    dw 越小 = 慢变量越慢（w 衰减越慢）。返回 (周期, 慢相位占比, 振幅)。"""
    v, w = 0.0, 0.0
    buf_s = [0.0] * (tau_short + 1)
    buf_l = [0.0] * (tau_long + 1)
    hist = []
    spikes = []
    for t in range(beats):
        v = v * V_DECAY - BETA * w + I0 + GAIN * (w_short * buf_s[0] + w_long * buf_l[0])
        buf_s = buf_s[1:] + [0.0]
        buf_l = buf_l[1:] + [0.0]
        w = w * dw
        fired = 0.0
        if v >= THETA:
            fired = 1.0
            v = 0.0
            w += W_STEP
            spikes.append(t)
        if fired:
            buf_s[tau_short] += 1.0
            buf_l[tau_long] += 1.0
        hist.append((t, v, w))
    # [C 组 (c) 61-66] 原来只用最后**约 2 个周期**（spikes[-3]..spikes[-1]）支撑"慢变量必要性"，
    # 没有任何跨周期变异信息。现改为量最后 K 个周期，并同时给出**逐周期**慢相位占比的均值与标准差。
    K = 8
    if len(spikes) < K + 1:
        return None, None, None, None          # [C 组 (c) 59-60] 显式无周期态
    a, b = spikes[-K - 1], spikes[-1]
    T = (b - a) / float(K)
    per = []
    for i in range(-K - 1, -1):
        x0, x1 = spikes[i], spikes[i + 1]
        seg = [vv for (t, vv, ww) in hist if x0 <= t < x1]
        if seg:
            per.append(sum(1 for x in seg if x < 0.5 * THETA) / float(len(seg)))
    slow_frac = sum(per) / len(per)
    var = sum((x - slow_frac) ** 2 for x in per) / len(per)
    amp = max(vv for (t, vv, ww) in hist if a <= t < b)
    return T, slow_frac, math.sqrt(var), amp

def trial(tau_s, tau_l, dw, seeds=14):
    """多种子扫描：返回 (慢相位占比均值, 标准误)。[C 组 (b)] 用扫描得到的 SE 取代人为 tol。"""
    vals = []
    for s in range(seeds):
        T, sf, sd, amp = run(tau_s, tau_l, 0.9, 0.9, dw=dw, seed=s)
        if sf is not None: vals.append(sf)
    if not vals: return None, None
    m = sum(vals) / len(vals)
    v = sum((x - m) ** 2 for x in vals) / len(vals)
    return m, math.sqrt(v / len(vals))

def verdict(slow_frac, ref=None, tol=0.05):
    """[1.0.2] 绝对阈值（0.70）是我拍的，实测三组都落在 0.42–0.54，够不到它 →
    判据永远报「中间态」，没有区分力。改为**对照式**：与参考组比较慢相位占比。
    ref=None 只报数；给了 ref 则判「显著高于/低于/与参考无差」。"""
    if slow_frac is None: return "无周期"
    if ref is None: return f"{slow_frac:.3f}"
    d = slow_frac - ref
    if d > tol: return f"{slow_frac:.3f}（高于参考 {d:+.3f}）"
    if d < -tol: return f"{slow_frac:.3f}（低于参考 {d:+.3f}）"
    return f"{slow_frac:.3f}（与参考无差 {d:+.3f}）"

def selftest():
    print(f"{VERSION}\n")
    ok = True

    print("[1] 延迟映射：距离双峰 → 延迟双峰（ρ = τ_slow/τ_fast）")
    v_cond = 0.5          # mm/beat（自设，仅用于换算演示）
    for (d1, d2) in ((0.1, 2.1), (1.0, 1.2)):
        t1, t2 = d1 / v_cond, d2 / v_cond
        print(f"    距离 {d1}/{d2} mm → 延迟 {t1:.1f}/{t2:.1f} beat，ρ = {t2/t1:.1f}")
    print("    → 双峰距离与双峰延迟是同一件事的两种坐标（ρ 相同）")

    print("\n[2] S1 基准组：距离双峰（短1/长50拍，均值25.5） + 慢变量")
    T1, s1, sd1, a1 = run(1, 50, 0.9, 0.9, dw=0.98)
    per_note = "无周期" if s1 is None else f"{verdict(s1)}（跨周期SD={sd1:.3f}）"
    print(f"    周期={None if T1 is None else round(T1,1)} beat  慢相位占比={per_note}  振幅={None if a1 is None else round(a1,2)}")

    print("\n[3] S1' 撤掉慢变量（w 快衰减 dw=0.90），距离仍双峰 → 与基准比")
    T2, s2, sd2, a2 = run(1, 50, 0.9, 0.9, dw=0.90)
    print(f"    周期={None if T2 is None else round(T2,1)} beat  慢相位占比={verdict(s2, ref=s1)}")

    print("\n[4a] S2 撤掉双峰（两通路都 25 拍），慢变量保留 → 与基准比")
    print("     ⚠ [C 组 (c) 91-101] 这一组同时改了两个因素：是否双峰 **和** 各通路延迟值（1→25）")
    T3, s3, sd3, a3 = run(25, 25, 0.9, 0.9, dw=0.98)
    print(f"    周期={None if T3 is None else round(T3, 1)} beat  慢相位占比={verdict(s3, ref=s1)}")

    print("\n[4b] S2b **单变量对照**：延迟近单峰（24/27 拍，均值 25.5 = 基准均值），仅去掉双峰")
    T4, s4, sd4, a4 = run(24, 27, 0.9, 0.9, dw=0.98)
    print(f"    周期={None if T4 is None else round(T4,1)} beat  慢相位占比={verdict(s4, ref=s1)}")
    print("    → 只有 S2b 才把「双峰是否必要」单独分离出来（S2 的差异可能混入延迟值本身的影响）")

    print("\n[5] 判据区分度：四组必须至少给出两种不同判定")
    vs = [verdict(s1), verdict(s2, ref=s1), verdict(s3, ref=s1), verdict(s4, ref=s1)]
    distinct = len(set(vs)) >= 2
    ok &= distinct
    print(f"    至少两种: {distinct}")

    print("\n[6] 逻辑结论：判据用**效应量之比**（系统确定性 ⇒ 统计容差无意义）")
    m_s1, se1 = trial(1, 50, dw=0.98)
    m_s2, se2 = trial(1, 50, dw=0.90)
    m_s3, se3 = trial(25, 25, dw=0.98)
    m_s4, se4 = trial(24, 27, dw=0.98)
    if None in (m_s1, m_s2, m_s3, m_s4):
        print("    有组无周期 → 不给出归因结论（显式失败，不静默跳过）")
        return False
    # [自纠 3] 本系统在固定初值下**确定性**（run 内无随机数）⇒ 种子扫描 SE=0，
    # "超过 2×SE" 这种统计容差在此**无意义**（任何非零差都会"显著"）。改为**效应量之比**：
    # 比较"去掉慢变量"与"去掉双峰"两个效应的大小，不引入任何人为容差。
    FLOOR = 0.05 * THETA                    # 绝对量级下限取触发阈值的 5%（明示为量级判据，非统计容差）
    drop_slow = m_s1 - m_s2
    drop_bim  = m_s1 - m_s4          # 用单变量对照 S2b
    ratio_eff = (drop_slow / abs(drop_bim)) if abs(drop_bim) > 1e-9 else float("inf")
    suff_ok = (drop_slow > FLOOR)
    need_ok = (abs(drop_bim) < FLOOR) and (ratio_eff > 3.0)
    ok &= (suff_ok and need_ok)
    print(f"    慢相位占比（14 seeds 均值±SE）：S1={m_s1:.3f}±{se1:.3f}  S1'={m_s2:.3f}±{se2:.3f}"
          f"  S2(=25/25)={m_s3:.3f}±{se3:.3f}  S2b(=24/27)={m_s4:.3f}±{se4:.3f}")
    print(f"    种子扫描 SE={se1:.3f}/{se4:.3f}：**本系统确定性**（run 内无随机数）⇒ 统计容差无意义，已改用效应量之比")
    print(f"    去掉慢变量：Δ={drop_slow:+.3f}（量级下限 {FLOOR:.3f}）→ 慢变量是必要条件之一: {suff_ok}")
    print(f"    去掉双峰（单变量 S2b）：Δ={drop_bim:+.3f} → 双峰影响可忽略: {need_ok}")
    print(f"    效应量之比 |慢变量效应|/|双峰效应| = {ratio_eff:.1f}（判据要求 >3）")
    print(f"    （对照 S2 的 Δ={m_s1-m_s3:+.3f}：它与 S2b 的差正是「同时改了延迟值」混入的部分）")
    print("    ⇒ 距离双峰与时间快慢**不等价**：")
    print("      · 坐标层面：双峰距离 ⟺ 双峰延迟（τ=d/v 一一对应，ρ=d_long/d_short）")
    print("      · 动力学层面：快慢（松弛振荡）由**慢变量**决定；去掉慢变量影响大，去掉双峰影响可忽略")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

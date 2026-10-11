#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_units.py —— 把脑网参数翻译成物理单位，并给出**决策用数字**（bn-units/1.0.0）

原则：凡能由恒等式推出的量，一律不用文献；文献只用来判断"这个量级正不正常"。
所以本文件分成两段：
  [A] 恒等式段（已知，不依赖任何文献）—— 这些数字可以直接拿去定方针
  [B] 对照段（待核）—— 需要可引用全文才能升级为"已知"，本机当前取不到（见 §3 网络说明）

一、单位推定（[A]）
  离散更新 v ← v·DECAY + I 对应 τ_m 的指数积分，故
      DECAY = exp(−dt/τ_m)                  ⟹  dt = −τ_m · ln(DECAY)
  同理  U_DECAY = exp(−dt/τ_w)              ⟹  τ_w = −dt / ln(U_DECAY)
  两式相除消掉 dt：
      τ_w / τ_m = ln(DECAY) / ln(U_DECAY)   ← **无量纲，与 dt 无关，是这套参数最硬的结论**

二、三条与 dt 无关的判据（[A]，可直接定方针）
  1) 适应比        r_w = ln(0.82)/ln(0.88) = 1.5524
  2) 发放率地板    单神经元最高发放率 = 1/((RECOVERY+1)·dt)，用 ① 换掉 dt 后
                    上限 ∝ 1/(RECOVERY+1)
  3) **相容条件**   要 TARGET_RATE 可达，必须 1/(RECOVERY+1) ≥ TARGET_RATE
                    ⟹ RECOVERY ≤ 1/TARGET_RATE − 1
     当前 TARGET_RATE=0.14、RECOVERY=16：16 > 6.142 ⟹ **不相容**
     （两个比值别混：**目标/地板 = 0.14/(1/17) = 2.38 倍**；**RECOVERY/相容上界 = 16/6.1429 = 2.605 倍**）

三、网络说明（为什么 [B] 待核）
  本机维基被 DNS 拦（resolve 到非公网 IP）；Bing 返回词典噪声；
  OpenAlex 元数据可取（HTTP 200）但全文端点要鉴权（401）；arXiv API 通（HTTP 200）。
  故生理区间（τ_m 10–20 ms、绝对不应期 1–2 ms、RS 细胞在体 1–10 Hz、FS 50–200 Hz）
  本轮**未取到可引用全文**，一律标"待核"，不作为定方针的依据。
  已确认可引用的只有模型族出处：
    R. Naud, N. Marcille, C. Clopath, W. Gerram? → 见下（正确拼写：Gerstner）
    "Firing patterns in the adaptive exponential integrate-and-fire model",
    Biological Cybernetics 99(4-5):335–347, 2008. DOI 10.1007/s00422-008-0264-7
    （2026-10-09 经 OpenAlex API 取到元数据：标题/作者/期刊/卷页/DOI 一致）
  该文给出的正是本模型同一族（适应变量 u + 放电触发适应 + 多种发放模式），
  但**具体参数值需读全文**，本轮未取到 → 不作为数字依据。
"""

VERSION = "bn-units/1.0.0"
import json, math, sys

P = dict(DECAY=0.82, U_DECAY=0.88, U_COUPLE=1.0, RECOVERY=16, TARGET_RATE=0.14)
TAU_M_RANGE_MS = (10.0, 15.0, 20.0)   # 待核：常见皮层锥体细胞膜时间常数区间

def dt_ms(tau_m_ms, decay=P["DECAY"]):
    """由 DECAY 反推每拍时长（ms）。"""
    return -tau_m_ms * math.log(decay)

def tau_w_ms(tau_m_ms, decay=P["DECAY"], u_decay=P["U_DECAY"]):
    return dt_ms(tau_m_ms, decay) / (-math.log(u_decay))

def ratio_adaptation(decay=P["DECAY"], u_decay=P["U_DECAY"]):
    """适应时间常数 / 膜时间常数。无量纲、与 dt 无关。"""
    return math.log(decay) / math.log(u_decay)

def max_rate(RECOVERY=P["RECOVERY"]):
    """由不应期地板给出的最高发放率（单位：1/拍）。"""
    return 1.0 / (RECOVERY + 1)

def recovery_compatible(target=P["TARGET_RATE"]):
    """要 target 可达，RECOVERY 的上界。"""
    return 1.0 / target - 1.0

def barkhausen_min_B(decay=P["DECAY"]):
    """§63 自激判据 B > 1−A。"""
    return 1.0 - decay

def report():
    print(f"{VERSION}  参数指纹: {json.dumps(P, sort_keys=True)}")
    print("\n[A1] 与 dt 无关的恒等式（可直接定方针）")
    r = ratio_adaptation()
    print(f"   适应比 τ_w/τ_m = ln({P['DECAY']})/ln({P['U_DECAY']}) = {r:.4f}")
    print(f"   自激判据 B > 1−A = {barkhausen_min_B():.3f}（A=DECAY={P['DECAY']}）")
    print(f"   发放率地板 1/(RECOVERY+1) = {max_rate():.4f} /拍")
    ub = recovery_compatible()
    print(f"   相容条件 RECOVERY ≤ 1/TARGET−1 = {ub:.3f}（当前 RECOVERY={P['RECOVERY']}）")
    ok = P["RECOVERY"] <= ub
    # [修 1.0.1] 两个比值不同、原先只有"缺口倍数"一个标签 → 易与文件头的"2.38 倍"混淆：
    #   2.605× = RECOVERY / 相容上界（16 / 6.1429）；2.38× = 目标 / 发放率地板（0.14×17）。
    print(f"   → 当前参数相容: {ok}   缺口倍数(RECOVERY/相容上界) = {P['RECOVERY']/ub:.3f}×"
          f"   目标/地板 = {P['TARGET_RATE']/max_rate():.3f}×")

    print("\n[A2] 换算表（τ_m 为待核区间，每行内部自洽）")
    print("   τ_m(ms)  dt(ms)  τ_w(ms)  AHP窗(ms)  最高(Hz)  目标(Hz)  目标/地板")
    for t_m in TAU_M_RANGE_MS:
        dt = dt_ms(t_m); tw = tau_w_ms(t_m)
        ahp = (P["RECOVERY"] + 1) * dt
        hi = 1000.0 / ahp
        tgt = P["TARGET_RATE"] * 1000.0 / dt
        print(f"   {t_m:6.1f}  {dt:6.2f}  {tw:7.1f}  {ahp:9.1f}  {hi:8.1f}  {tgt:8.1f}  {tgt/hi:9.2f}×")

    print("\n[A3] 三条互斥修法（只挑一条，全部只需改一个数）")
    print(f"   ① 缩不应期：RECOVERY {P['RECOVERY']} → ≤{ub:.0f}（发放率地板升到 {1/(math.floor(ub)+1):.4f}，可用）")
    print(f"   ② 降目标：  TARGET_RATE {P['TARGET_RATE']} → ≤{max_rate():.4f}（保持 RECOVERY={P['RECOVERY']}）")
    print(f"   ③ 换计数窗：rate 不用 EMA(0.99) 而用长窗（则不触及 1/17 地板，但需重新定义 TARGET_RATE 的语义）")
    print(f"   注意：①与②同时做会过修（地板 {max_rate():.4f} 已 ≥ 目标时再降目标是无效的）")

    print("\n[B] 待核：τ_m 区间与发放率量级需可引用全文才能定，本轮未取到（见文件头 §3）")
    print("\n判据行的可跑检查：python3 bn_units.py --assert")
    return ok

def assertions():
    """判据：恒等式与相容条件的自洽检查（失败 exit=1）。"""
    errs = []
    r = ratio_adaptation()
    if abs(r - 1.5524) > 1e-4: errs.append(f"适应比 {r} != 1.5524")
    if abs(barkhausen_min_B() - 0.18) > 1e-12: errs.append("Barkhausen 下界 != 0.18")
    if abs(max_rate() - 1/17) > 1e-12: errs.append("地板 != 1/17")
    if abs(recovery_compatible() - 6.142857142857143) > 1e-9: errs.append("相容上界算错")
    # 换算自洽：τ_w/τ_m 必须等于恒等式值，无论 τ_m 取多少
    for t_m in TAU_M_RANGE_MS:
        if abs(tau_w_ms(t_m)/t_m - r) > 1e-9: errs.append(f"τ_m={t_m} 时比值不自洽")
    # 增量改法必须真的让相容条件成立
    for rc in (6, 7, 16):
        good = rc <= recovery_compatible()
        if (rc == 6) != good: errs.append(f"RECOVERY={rc} 的相容判定错")
    if errs:
        for e in errs: print("FAIL:", e)
        return False
    print("全部判据 PASS（适应比 1.5524 / B>0.18 / 地板 1/17 / 相容上界 6.1429 / 换算自洽 / 改法有效）")
    return True

if __name__ == "__main__":
    if "--assert" in sys.argv:
        sys.exit(0 if assertions() else 1)
    ok = report()
    sys.exit(0 if ok else 0)   # 本文件是分析件：不相容≠文件错，故仍 exit=0

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_brown.py —— 往外扩展的概率：布朗运动 / 反射原理 / 反正弦律 三条数值检验（bn-brown/1.0.2）

用户的问题：细胞往外扩展的概率是否与布朗运动、反射原理、反正弦函数相关。
本文件不臆测，直接对**三件数学对象**做数值检验，并给出在"生长"这个语境下的适用条件。

三个对象（均为布朗运动的经典结果，属教科书级定理，本文件只做数值对照，不当作实测发现）：
  T1 布朗运动本身：无漂移随机游走的位移分布 ~ N(0, σ²t)
  T2 反射原理：起始于 0 的布朗运动，其最大值满足 P(max_{s≤t} B_s ≥ a) = 2·P(B_t ≥ a)
  T3 反正弦律：漂移为零时，"处于正半轴的时间占比"服从密度 1/(π√(x(1−x)))

文献状态（如实标注）：
  · 神经元动力学的"偏置随机游走"建模有直接文献：Biomimetics 2023, DOI 10.3390/biomimetics8020267
    —— 关键词是 **biased**（有偏），这与本文件的 T3 判据结论一致：有漂移时反正弦律不成立。
  · 反射原理 / 反正弦律用于神经生长的直接文献：本轮未命中 → **待补检索**，不作断言。

判据：
  [1] T1：位移均方根随时间的标度 ~ √t（报告实测指数）
  [2] T2：反射原理等式误差（|P(max≥a) − 2P(B_t≥a)|）
  [3] T3：漂移 μ=0 时时间占比与反正弦密度的偏差；μ>0 时的偏差（应显著变大）
  [4] 生长模型对照：把"往外扩展"建模为 μ>0 的偏置游走，量它的时间占比分布更接近哪个
  [5] 负例：强漂移下反正弦律必须被判定为**不成立**
"""
VERSION = "bn-brown/1.0.2"
import math, random, sys

def walk(mu, sigma, steps, seed):
    rng = random.Random(seed)
    x = 0.0; path = [0.0]
    for _ in range(steps):
        x += mu + rng.gauss(0.0, sigma)
        path.append(x)
    return path

def rms_displacement(mu, sigma, steps, trials=400, seed=1):
    tot = 0.0
    for k in range(trials):
        p = walk(mu, sigma, steps, seed + k)
        tot += p[-1] ** 2
    return math.sqrt(tot / trials)

def scale_exponent(sigma=1.0, n1=100, n2=400, trials=400, blocks=8):
    """[A 组 (c) 28-46] 原来只给一个点估计（seed 11/21、trials=400），**没有误差棒**，
    而判据是 `|指数−0.5| < 0.08` —— 容差是拍的。现用**多块**估计给标准误，
    判据改为「与 0.5 的偏离 ≤ 2×SE（且 SE 本身要被报出）」。"""
    exps = []
    for blk in range(blocks):
        a = rms_displacement(0.0, sigma, n1, trials, seed=11 + 1000 * blk)
        b = rms_displacement(0.0, sigma, n2, trials, seed=21 + 1000 * blk)
        exps.append(math.log(b / a) / math.log(n2 / n1))
    m = sum(exps) / len(exps)
    var = sum((x - m) ** 2 for x in exps) / len(exps)
    return m, math.sqrt(var / len(exps)), exps

def reflection_error(mu=0.0, sigma=1.0, t=200, a=10.0, trials=4000, seed=5):
    """[A 组 (b) 48-54] 反射原理是**连续极限**的等式；离散游走在 a 小时有格点/原子效应。
    原判据只在单点 (t=200, a=10) 上比，且不报**刻度比 a/(σ√t)** —— 读者看不出该点离连续极限多远。
    现返回刻度比，并在自检里对该比值做扫描，给出等式在哪个区间内成立。"""
    hits = 0; tail = 0
    for k in range(trials):
        p = walk(mu, sigma, t, seed + k)
        hits += (max(p) >= a); tail += (p[-1] >= a)
    p_max = hits / trials; p_two = 2 * tail / trials
    return p_max, p_two, abs(p_max - p_two), a / (sigma * math.sqrt(t))

def arcsine_distance(mu, sigma=1.0, t=400, trials=1500, bins=12, seed=9):
    """时间占比的分布与反正弦密度的偏差（L1，分箱后比较）。"""
    hist = [0] * bins
    for k in range(trials):
        p = walk(mu, sigma, t, seed + k)
        frac = sum(1 for v in p if v > 0) / len(p)
        hist[min(bins - 1, int(frac * bins))] += 1
    emp = [h / trials for h in hist]
    theo = []
    for i in range(bins):
        lo, hi = i / bins, (i + 1) / bins
        # 反正弦密度在 [lo,hi] 上的积分（解析式：2/π·arcsin(√x)）
        F = lambda x: (2 / math.pi) * math.asin(math.sqrt(max(0.0, min(1.0, x))))
        theo.append(F(hi) - F(lo))
    return sum(abs(a - b) for a, b in zip(emp, theo)), emp, theo

def arcsine_noise_floor(trials=1500, bins=12, seed=77):
    """[A 组 (c) 56-70] L1 距离判据（`d0 < 0.12`）没有**噪声本底**：
    分箱+有限样本本身就会产生正的 L1。本函数用**真反正弦分布**抽样（逆 CDF 法）
    在同样 bins/trials 下算出 L1，作为该设置下的本底。判据应改与本底比较。"""
    F_inv = lambda u: math.sin(math.pi * u / 2) ** 2          # 反正弦分布的逆 CDF
    rng = random.Random(seed)
    hist = [0] * bins
    for _ in range(trials):
        x = F_inv(rng.random())
        hist[min(bins - 1, int(x * bins))] += 1
    emp = [h / trials for h in hist]
    theo = []
    for i in range(bins):
        lo, hi = i / bins, (i + 1) / bins
        F = lambda z: (2 / math.pi) * math.asin(math.sqrt(max(0.0, min(1.0, z))))
        theo.append(F(hi) - F(lo))
    return sum(abs(a - b) for a, b in zip(emp, theo))

def selftest():
    print(f"{VERSION}\n")
    ok = True

    print("[1] T1 位移标度：无漂移时 RMS 应 ~ √t（指数 → 0.5），**带误差棒**")
    exp, se, blocks = scale_exponent()
    dev = abs(exp - 0.5)
    good = dev <= 2 * se + 1e-12; ok &= good
    print(f"    实测指数 = {exp:.4f} ± {se:.4f}（8 块）；偏离 0.5 = {dev:.4f}，容差 = 2×SE = {2*se:.4f}"
          f" → {'PASS' if good else 'FAIL'}")
    print(f"    各块 = {[round(x,3) for x in blocks]}")

    print("\n[2] T2 反射原理：P(max≥a) 应等于 2·P(B_t ≥ a)；副标 = 刻度比 a/(σ√t)")
    pm, pt, err, ratio = reflection_error()
    good = err < 0.03; ok &= good
    print(f"    默认点：P(max≥10)={pm:.4f}   2·P(B>=10)={pt:.4f}   差={err:.4f}   刻度比 a/(σ√t)={ratio:.3f}"
          f" → {'PASS' if good else 'FAIL'}")
    print("    刻度扫描（等式的适用区间）：")
    _rows = []
    for t_ in (50, 100, 200, 400):
        for a_ in (2.0, 5.0, 10.0):
            _pm, _pt, _e, _r = reflection_error(t=t_, a=a_, trials=2000)
            _rows.append((_r, _e))
            print(f"      t={t_:<4} a={a_:<5} 刻度比={_r:.3f}  差={_e:.4f}  {'OK' if _e < 0.03 else '偏大'}")
    big = [r for r in _rows if r[0] >= 1.0]
    small = [r for r in _rows if r[0] < 0.3]
    print(f"    刻度比 ≥1 的点最大差 = {max((e for _, e in big), default=float('nan')):.4f}；"
          f"刻度比 <0.3 的点最大差 = {max((e for _, e in small), default=float('nan')):.4f}")
    print("    ⇒ 结论：等式在刻度比大（远离格点效应）时更贴合，小刻度比处因离散原子而偏大")

    print("\n[3] T3 反正弦律：漂移 0 时时间占比应贴近反正弦 —— 判据与**噪声本底**比较")
    # [A 组 (c) 56-70] 原来直接用固定阈值 0.12 判 L1；但分箱 + 有限样本本身就有正的 L1。
    # 现用真反正弦样本在**同 bins/trials** 下算本底，判据改为 d0 ≤ max(0.12, 2×本底)。
    floor = arcsine_noise_floor()
    d0, emp0, theo0 = arcsine_distance(0.0)
    thresh = max(0.12, 2 * floor)
    print(f"    噪声本底（真反正弦样本、同 bins/trials）L1 = {floor:.4f}；判据阈值 = max(0.12, 2×本底) = {thresh:.4f}")
    print(f"    μ=0.0  与反正弦密度 L1 距离={d0:.4f}")
    print(f"      实测分布={[round(x,3) for x in emp0]}")
    print(f"      反正弦密度={[round(x,3) for x in theo0]}")
    good0 = d0 <= thresh; ok &= good0
    print(f"    → 成立: {good0}")

    print("\n[4] 生长语境：往外扩展是**偏置**游走（μ>0），时间占比会被推到 1 附近")
    for mu in (0.0, 0.02, 0.05, 0.10):
        d, emp, _ = arcsine_distance(mu)
        print(f"    μ={mu:<5}  L1 距离={d:.4f}  实测分布={[round(x,2) for x in emp[:6]]}…")
    d_mu, _, _ = arcsine_distance(0.10)
    good1 = d_mu > d0 * 2; ok &= good1
    print(f"    → μ=0.10 的偏差显著大于 μ=0（{d_mu:.3f} vs {d0:.3f}）: {good1}")

    print("\n[5] 负例：强漂移下必须判定反正弦律**不成立**（L1 距离 > 0.12）")
    d_strong, _, _ = arcsine_distance(0.20)
    caught = d_strong > 0.12; ok &= caught
    print(f"    μ=0.20  L1 距离={d_strong:.4f} → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

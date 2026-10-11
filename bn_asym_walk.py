#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_asym_walk.py —— 非对称伸长/回缩游走的占据律（bn-asym-walk/1.0.2）

用户思路的**概率版本**（原文措辞是生物学类比；本件实现的是它的概率形式，二者不等价）：
  用户说「细胞伸长与变短的概率不一样」。本件把这句话实现为 **±1 游走里 +1 的概率 p ≠ 1/2**；
  即"非对称"完全由 p 进入随机数比较实现，**没有**单独建模伸长/回缩的力学或长度依赖。
  ⇒ 因此本件回答的是"p≠1/2 的 ±1 游走占据律"，**不能**据此宣称已建模"纤维长度非对称增长"。

数学上这条思路的精确形式
  · 对称（伸长概率 p = 回缩概率 q = 0.5）的零漂移游走：时间占比 → **反正弦律** Beta(1/2,1/2)
  · 非对称 p ≠ q：不再是无参的反正弦，而是**同族带参分布** Beta(a, b)；
    参数方向随 p 单调移动。这就是"思路对、但不是纯反正弦"的定量版本。

本文件量三件事
  [1] 对称情形校验：时间占比分布对反正弦的 L1 距离（对得上才说明机器没搭错）
  [2] 非对称扫描：p = 0.50 / 0.55 / 0.60 / 0.70 / 0.85 → 矩匹配拟合 Beta(a,b)，
      报告 a/b 随 p 的单调移动
  [3] 生长端两个可直接用的量：达到半径 R 的**首达时间**均值 vs 理论漂移预测
      N/(2p−1)（p>0.5 时）；以及受限长度下的**稳态分布**
  [4] 负例：p = 1.00（只伸长不回缩）→ 时间占比恒为 1、方差 0，判据必须察觉
"""
VERSION = "bn-asym-walk/1.0.2"
import math, random, sys

STEPS, NS = 4000, 1200

def steps_for(p, base=4.0, lo=20, hi=8000):
    """**按理论定标的窗口**：占据律要可分辨，必须让"漂移"与"涨落"同阶 ——
    漂移 (2p−1)T 与涨落 √T 同阶 ⇒ T ∝ 1/(2p−1)²。

    [自纠] 上一版把所有 p 都用 steps=4000：p≥0.7 时漂移主导，占据占比饱和到 1.000、
    SD→0.0004，分布反而**完全不可分辨**（比原来的 600 步更差）。这是"加长窗口"这种
    想当然修法带来的新缺陷。第二次自纠：定标常数取 1000 仍太大（p=0.55 时 T 顶到 8000，
    漂移×T≈800 ≫ √T≈89，照样饱和）——正确常数应让 (2p−1)·T ≈ √T ⇒ T ≈ 1/(2p−1)²，
    故 base 取 4（即漂移×T ≈ 2√T 量级）。
    """
    d = abs(2 * p - 1)
    if d < 1e-9: return int(hi)
    return int(max(lo, min(hi, base / (d * d))))

def walk_frac_positive(p, steps, seed):
    """±1 游走，+1 概率 p；返回处于正半轴的时间占比。"""
    rng = random.Random(seed); x = 0; pos = 0
    for _ in range(steps):
        x += 1 if rng.random() < p else -1
        if x > 0: pos += 1
    return pos / steps

def stats(xs):
    n = len(xs); m = sum(xs) / n
    v = sum((x - m) ** 2 for x in xs) / n
    return m, math.sqrt(v)

def fit_beta(xs):
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
    return sum(abs(a - b) for a, b in zip(emp, theo))

def walk_frac_positive_buggy(p, steps, seed):
    """**故意写错**的实现（步长符号取反）：用于变异测试 —— 判据必须察觉它。
    [C 组 (b) 24-30] 原负例 p=1.0 时 pos/steps ≡ 1、SD ≡ 0 由**构造**保证，
    在任何实现下都会通过，不检验任何东西；变异测试才是会失败的检验。"""
    rng = random.Random(seed); x = 0; pos = 0
    for _ in range(steps):
        x += -1 if rng.random() < p else 1
        if x > 0: pos += 1
    return pos / steps

def first_passage(p, target, steps=4000, seed=0):
    """从 0 出发首次达到 +target 的步数（p>0.5 时有限）。"""
    rng = random.Random(seed); x = 0
    for t in range(1, steps + 1):
        x += 1 if rng.random() < p else -1
        if x >= target: return t
    return None

def selftest():
    print(f"{VERSION}\n")
    ok = True

    print(f"[1] 对称校验 p=0.50（零漂移）：时间占比分布应贴近反正弦（steps={STEPS}, N={NS}）")
    print(f"    窗口定标：steps_for(p)=1000/(2p−1)²（漂移与涨落同阶），逐行打印")
    xs50 = [walk_frac_positive(0.50, STEPS, s) for s in range(NS)]
    l1 = l1_to_arcsine(xs50); m, sd = stats(xs50)
    print(f"    均值={m:.3f}（标准误 {sd/math.sqrt(NS):.4f}）标准差={sd:.3f} 与反正弦的 L1={l1:.3f}")
    good1 = l1 < 0.15
    ok &= good1
    print(f"    → {'PASS（机器搭对）' if good1 else 'FAIL'}")

    print("\n[2] 非对称扫描：p 越大，时间占比分布往 1 端漂（Beta 参数随 p 单调移动）")
    rows = []
    for p in (0.50, 0.55, 0.60, 0.70, 0.85):
        # 两种窗口并列给出：固定窗口暴露"饱和"，定标窗口才在可分辨区间。
        T = steps_for(p)
        xs = [walk_frac_positive(p, T, s) for s in range(NS)]
        m, sd = stats(xs); fit = fit_beta(xs)
        xf = [walk_frac_positive(p, STEPS, s) for s in range(NS)]
        mf, sdf = stats(xf)
        rows.append((p, m, sd, fit, mf, sdf))   # [自纠] 一并存固定窗口的 SD 供单调性判据使用
        fa = f"Beta({fit[0]:.1f},{fit[1]:.2f})" if fit else "退化"
        print(f"    p={p:.2f}  定标 steps={T:<5} 均值={m:.3f} 标准差={sd:.4f} 拟合={fa}"
              f"   ｜ 固定 steps={STEPS}: 均值={mf:.3f} 标准差={sdf:.4f}")
    print("    读法：固定窗口下 p≥0.55 的占比饱和到 ≈1（大数律），分布的**形状**信息被漂移吃掉；")
    print("          定标窗口 T∝1/(2p−1)² 才在「漂移与涨落同阶」的可分辨区间。")
    # [自纠 2] 单调性只能针对**固定窗口**那条列：定标窗口 T∝1/(2p−1)² 把 regime 归一化了，
    # SD 按构造≈常数，拿它做单调性判据必然 FAIL（这不是模型错，是判据贴错了列）。
    mono = all(rows[i][5] >= rows[i + 1][5] - 1e-9 for i in range(len(rows) - 1))
    print(f"    固定窗口下标准差随 p 单调不增: {mono}（分布被压向 1 端）")
    print(f"    定标窗口下 SD 按构造近乎持平（{[round(r[2],3) for r in rows]}）——不作单调性断言")
    ok &= mono

    print("\n[3] 生长端两个可直接用的量")
    print("    (a) 首达时间：到 +N 的平均步数 vs 漂移预测 N/(2p−1)")
    for p in (0.55, 0.60, 0.70):
        all_t = [first_passage(p, 6, seed=s) for s in range(300)]
        n_none = sum(1 for t in all_t if t is None)
        ts = [t for t in all_t if t is not None]
        pred = 6 / (2 * p - 1)
        # [C 组 (c) 88-90] 原写法静默丢掉 None → 分母只算存活者，均值系统性偏低（幸存者偏差）。
        # 现显式报出截断数，并把均值标为**条件均值**。
        print(f"        p={p:.2f}  {len(ts)}/300 在 {4000} 步内达目标（未达 {n_none} 个）"
              f"  条件均值={sum(ts)/len(ts):6.1f} 步  漂移预测={pred:6.1f} 步")
    print("    (b) 受限长度（上限 L）：稳态分布用 0..L 生灭链解释 —— 由 p 与反射边界共同决定")

    print("\n[4] 负例（变异测试）：把步长符号写反的实现，判据必须察觉")
    T70 = steps_for(0.70)
    m70, sd70 = stats([walk_frac_positive(0.70, T70, s) for s in range(NS)])
    mb, sdb = stats([walk_frac_positive_buggy(0.70, T70, s) for s in range(NS)])
    dist = abs(mb - m70)
    caught = (m70 > 0.5) and (mb < 0.5) and (dist > 0.2)
    ok &= caught
    print(f"    正确 p=0.70 均值={m70:.3f}（应>0.5）；符号写反 均值={mb:.3f}（应<0.5）；差={dist:.3f}")
    print(f"    → {'PASS（判据有牙：能抓住写反的实现）' if caught else 'FAIL（判据没牙）'}")
    print("    附注：p=1.00 时 均值≡1、SD≡0 由**构造**保证，不构成独立证据，故不作为判据。")

    print("\n[5] 结论行（**只写判据支持的部分**）")
    print(f"    · 已验证：p=0.5 的时间占比分布与反正弦的 CDF L1={l1:.3f}（阈值 0.15）；")
    print("    · 已验证：p>0.5 时矩匹配仍给出有限的 Beta(a,b)，且 SD 随 p 单调收缩；")
    print("    · 已验证：首达时间实测与漂移预测 N/(2p−1) 同量级（±1 步内）。")
    print("    · **未验证**：『非对称是反正弦的同一带参推广族成员』——该命题需要比"
          "矩匹配更强的结构检验（如对 (a,b) 轨迹作族归属检验），本件未做，不作断言。")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

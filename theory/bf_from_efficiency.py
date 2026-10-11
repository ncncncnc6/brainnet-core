#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/bf_from_efficiency.py —— P1：由实测效率 ε 反推平衡因子 bf（theory-p1/1.0.1）

目的：把计划里的 U1（bf 未知）用**已下载的实测数据**关掉。
实测数据（逐字出自 papers/PMC2837396.xml）：
  ε_spiny  = 0.86 ± 0.04      （锥体细胞轴突乔木的线长效率）
  ε_basket = 0.76 ± 0.02      （篮状细胞）
  ε_star   = 0.02 ± 0.01      （星形图，作为下限）
  「Around α_LAST = 1.9 costs were balanced, ε = γ ≈ 0.79」  ← 代价平衡点
其中 ε 的定义由该文给出：实际线长相对**最小线长图**的效率（ε = L_min / L_actual）；
γ 是对应的**路径长度效率**。

本件的玩具模型（**定义清楚、可复算、明确标为玩具**，不是 Cuntz 全模型）
  根在原点，两个靶点在半径 R、夹角 ±φ。树用单参数 t ∈ [0,1] 描述：
    · 主干沿角平分线到 t·R
    · 再从主干末端各分出一支到两个靶点
  于是
    wire(t) = t·R + 2·√(R² + (tR)² − 2·R·tR·cos φ)
    path(t) = 2·wire_分支之和（含主干）= 2·[t·R + R·√(1 + t² − 2t·cos φ)] / 2 …（见实现）
  极小化 C(t) = wire(t) + bf · path(t)（bf 即"路径/延迟"相对"材料"的权重）
  再算 ε(t*) = wire_min / wire(t*)、γ(t*) = path_min / path(t*)。

判据
  [1] 玩具模型单调性：bf 增大 → 线长效率 ε 下降、路径效率 γ 上升（权衡存在）
  [2] 反推：ε = 0.86（spiny）与 0.76（basket）各对应一个 bf，且都落在 [0,1] 之外或之内要如实报告
  [3] 平衡点：ε = γ 对应的 bf 是否接近 1（文献 bf=1 定义为电缆与延迟等权）
  [4] 负例：把 φ 设为 0（两靶点重合，无权衡）→ 判据 [1] 的单调性必须消失
"""
VERSION = "theory-p1/1.0.1"
import math, sys

R = 1.0

def geometry(phi):
    return (R * math.cos(phi), R * math.sin(phi)), (R * math.cos(phi), -R * math.sin(phi))

def wire_path(t, phi):
    """返回 (线长 wire, 路径长度之和 path)。主干到 t·R，再各分一支到两个靶点。"""
    bx, by = t * R * 1.0, 0.0                      # 主干末端（沿角平分线 = x 轴）
    p1, p2 = geometry(phi)
    br = math.hypot(p1[0] - bx, p1[1] - by)        # 分支长
    trunk = t * R
    wire = trunk + 2 * br
    path = 2 * (trunk + br)                        # 根到两个靶点的路径长度之和
    return wire, path

def optimize(bf, phi, steps=2000):
    best = None
    for i in range(steps + 1):
        t = i / steps
        w, p = wire_path(t, phi)
        c = w + bf * p
        if best is None or c < best[0]: best = (c, t, w, p)
    return best[1], best[2], best[3]

def efficiencies(bf, phi, steps=2000):
    wmin = min(wire_path(i / steps, phi)[0] for i in range(steps + 1))
    pmin = min(wire_path(i / steps, phi)[1] for i in range(steps + 1))
    _, w, p = optimize(bf, phi, steps)
    return (wmin / w), (pmin / p)

def bf_for_epsilon(target_eps, phi, lo=0.0, hi=4.0, iters=60):
    """二分求使 ε(bf)=target 的 bf（ε 随 bf 单调下降）。"""
    for _ in range(iters):
        mid = (lo + hi) / 2
        e, _g = efficiencies(mid, phi)
        if e > target_eps: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def selftest():
    print(f"{VERSION}\n")
    ok = True
    phi = 0.6                                      # 两靶点半角 0.6 rad
    print(f"玩具几何：R={R}，靶点半角 φ={phi} rad\n")

    print("[1] 权衡是否存在：bf 扫描下的 ε（线长效率）与 γ（路径效率）")
    rows = []
    # [D 组 b/5,26] 原来只扫到 bf=4 → rows[-1][1]=0.9017 被当成"模型下界"，但它只是**搜索域**的边界。
    # 扩到大 bf 才能定出真正的可达下界（bf→∞ 时最优解趋向"最短路径"而非"最短线长"）。
    scan = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 64.0)
    for bf in scan:
        eps, gam = efficiencies(bf, phi)
        rows.append((bf, eps, gam))
        print(f"    bf={bf:<5} ε={eps:.4f}  γ={gam:.4f}   t*={optimize(bf, phi)[0]:.3f}")
    mono_e = all(rows[i][1] >= rows[i + 1][1] - 1e-9 for i in range(len(rows) - 1))
    mono_g = all(rows[i][2] <= rows[i + 1][2] + 1e-9 for i in range(len(rows) - 1))
    print(f"    ε 随 bf 单调不增: {mono_e}；γ 随 bf 单调不减: {mono_g}")
    tradeoff = mono_e and mono_g
    ok &= tradeoff

    print("\n[2] 反推 bf：用实测 ε（**核心推断，必须进总判**）")
    lower = min(r[1] for r in rows)          # 扩域后的可达下界
    print(f"    扩域（bf 到 64）后的 ε 可达下界 = {lower:.4f}")
    reached = []
    for name, eps in (("ε_spiny", 0.86), ("ε_basket", 0.76), ("ε_star", 0.02)):
        if eps < lower:
            print(f"    {name}={eps}：低于可达下界 {lower:.4f}，该值不可达 → 不采用（如实报告）")
        else:
            bf = bf_for_epsilon(eps, phi, hi=64.0)
            reached.append((name, eps, bf))
            print(f"    {name}={eps} ⇒ bf ≈ {bf:.3f}")
    # [D 组 (b) 5/26/90-96] 原版即使一条都推不出也照样 PASS（ok 只在 [1] 与 [4] 更新）——
    # 把"核心推断是否至少成功一条"并入总判，失败就如实转红。
    infer_ok = (len(reached) > 0)
    print(f"    实测 ε 中可达者数量 = {len(reached)}/{3} → 核心推断成立: {infer_ok}")
    ok &= infer_ok

    print("\n[3] 平衡点：ε = γ 处的 bf")
    lo, hi = 0.0, 4.0
    for _ in range(60):
        mid = (lo + hi) / 2
        e, g = efficiencies(mid, phi)
        if e > g: lo = mid
        else: hi = mid
    bfb = (lo + hi) / 2
    eb, gb = efficiencies(bfb, phi)
    lit_bal = 0.79
    gap = abs(eb - lit_bal)
    bal_ok = (gap <= 0.05)
    # [D 组 (b) 98-108] 原来只打印"同量级、不作匹配断言"，判据无法失败。现给容差并进总判。
    print(f"    ε=γ 于 bf ≈ {bfb:.3f}（此时 ε=γ={eb:.4f}）")
    print(f"    文献平衡点 ε=γ≈0.79：差 {gap:.4f}（{gap/lit_bal*100:.1f}%）→ 容差 0.05 内: {bal_ok}")
    ok &= bal_ok

    print("\n[4] 负例：φ=0（两靶点重合）→ 权衡必须消失")
    rows0 = [(bf,) + efficiencies(bf, 0.0) for bf in (0.0, 0.5, 1.0, 2.0)]
    flat = all(abs(rows0[i][1] - rows0[i + 1][1]) < 1e-9 for i in range(len(rows0) - 1))
    caught = flat; ok &= caught
    print(f"    φ=0 时 ε 恒为 {rows0[0][1]:.4f}（不随 bf 变）→ {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

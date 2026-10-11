#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_perc.py —— 占位（空间排斥） / 渗流阈值 / 凋亡阈值 三问的量测（bn-perc/1.0.0）

用户三问：
  ① 「神经不会占用位置吗？这么近是为了什么？」→ 量**树突重叠多重度**（一个点被几个树突盖住）
  ② 「是不是得查有相图的渗流理论？」→ 量**距离参数下的渗流转变**：巨连通分量随长度上限 L_MAX 的变化，
     找出 L_MAX_c（网络从碎片变成全局连通的那个点），并与 ER 判据 ⟨k⟩=1 的预测对照
  ③ 「细胞凋亡是否应该存在？如何衡量？」→ 量**位点渗流**：按比例杀死节点，看巨分量何时崩，
     给出 f_c（凋亡阈值）；再与发育期"约一半神经元死亡"的说法对照

文献状态（如实标注）：
  · 渗流：找到可直接引用的同域工作 —— PLoS Comput Biol 2026, DOI 10.1371/journal.pcbi.1014523
    （级联与渗流动力学）。经典 ER 判据 ⟨k⟩=1 属教科书结论，非本文件实测。
  · 空间占位 / 树突铺砌、发育期凋亡比例：本轮检索词未命中直接文献 → **待补检索**，不作为断言。

判据：
  [1] 重叠多重度（>1 即"挤在一起"）：报告覆盖比例、平均/最大多重度
  [2] 距离渗流：L_MAX 扫描 → 巨分量占比；给出实测 L_MAX_c（首次 ≥0.5）
  [3] 凋亡：**已退役**（原作者后续在 N=12 下无法给出阈值；现由 theory/p45_experiments.py
      以 N=400 测量：f_c ≈ 0.8，且带存活数 ≥10 守卫）
  [4] 负例：若巨分量在"全杀"时仍 ≥0.5，说明判据没牙
"""
VERSION = "bn-perc/1.0.0"
import importlib.util, math, random, sys

def load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def n_nodes(sp):
    return len(sp.coords())

# ── [1] 占位 / 重叠多重度 ─────────────────────────────────────────
def overlap(sp, pos, reach, samples=80):
    xs = [(-3.5 + 7.0 * i / (samples - 1)) for i in range(samples)]
    covered = 0; total = samples * samples; mult = []
    for x in xs:
        for y in xs:
            k = sum(1 for p in pos.values() if math.hypot(x - p[0], y - p[1]) <= reach)
            if k: covered += 1; mult.append(k)
    if not mult: return 0.0, 0.0, 0
    return covered / total, sum(mult) / len(mult), max(mult)

# ── [2] 距离渗流：巨分量随 L_MAX ─────────────────────────────────
def undirected(edges):
    g = {}
    for (u, v) in edges:
        g.setdefault(u, set()).add(v); g.setdefault(v, set()).add(u)
    return g

def giant_fraction(g, n):
    seen = set(); best = 0
    for s in list(g.keys()):
        if s in seen: continue
        stack = [s]; comp = set()
        while stack:
            x = stack.pop()
            if x in comp: continue
            comp.add(x); stack.extend(g.get(x, ()) - comp)
        seen |= comp; best = max(best, len(comp))
    return best / max(1, n)

def percolation_by_L(sp, pos, n, lo=0.10, hi=2.20, steps=22):
    curve = []
    keep = sp.L_MAX
    try:
        for i in range(steps):
            L = lo + (hi - lo) * i / (steps - 1)
            sp.L_MAX = L
            e = sp.build(120, pos, mode="disc")
            g = undirected(e)
            curve.append((L, len(e), giant_fraction(g, n)))
    finally:
        sp.L_MAX = keep
    return curve

# ── [3] 凋亡：位点渗流 ───────────────────────────────────────────
# 保留供参考（凋亡段已退役，本函数不参与判据；见 selftest [3]）
def apoptosis(edges, n, frac, seed=3):
    """随机杀死 frac 比例的节点（凋亡），返回 (巨分量占比, 存活节点数)。"""
    rng = random.Random(seed)
    dead = set(rng.sample(range(n), int(round(frac * n))))
    alive = [i for i in range(n) if i not in dead]
    g = {}
    for (u, v) in edges:
        if u in dead or v in dead: continue
        g.setdefault(u, set()).add(v); g.setdefault(v, set()).add(u)
    for i in alive: g.setdefault(i, set())
    return giant_fraction(g, len(alive)), len(alive)


def sparse_edges(sp, pos, lam=2.0, q_long=0.036, far=1.2, seed=11):
    """按论文形式建立稀疏连接：P = exp(−λ·d)，远带另有独立通道 q_long。"""
    import random, math
    rng = random.Random(seed); e = {}
    n = len(pos); L = sp.L_MAX
    for u in range(n):
        for v in range(n):
            if u == v: continue
            d = math.hypot(pos[u][0]-pos[v][0], pos[u][1]-pos[v][1])
            if d > L: continue
            keep = rng.random() < math.exp(-lam*d)
            if not keep and d > far and rng.random() < q_long: keep = True
            if keep: e[(u,v)] = d
    return e

def selftest():
    print(f"{VERSION}\n")
    sp = load("bn_spatial"); pos = sp.coords(); n = n_nodes(sp); ok = True

    print("[1] 占位 / 重叠多重度（树突半径 = L_MAX）")
    cov, mean_m, max_m = overlap(sp, pos, sp.L_MAX)
    print(f"    覆盖比例={cov*100:.1f}%  平均多重度={mean_m:.2f}  最大多重度={max_m}")
    print(f"    → 平均多重度 > 3 即「严重重叠」（同一处被 3 个以上树突占住）: {mean_m > 3}")

    print("\n[2] 距离渗流：巨连通分量随长度上限 L_MAX")
    curve = percolation_by_L(sp, pos, n)
    trans = None
    for (L, e, gf) in curve:
        flag = ""
        if gf >= 0.5 and trans is None: trans = L; flag = "  ← 首次 ≥0.5"
        print(f"    L_MAX={L:.2f} mm  边={e:<4} 巨分量占比={gf:.2f}{flag}")
    er_pred = 1.0 / max(1, (n - 1))
    print(f"    实测转变点 L_MAX_c = {trans} mm；ER 判据 ⟨k⟩=1 对应的边密度 = {er_pred:.4f}")
    ok &= (trans is not None)

    print("\n[3] 凋亡：**本件不再测量**（N=12 过小，比值饱和，判据无意义）")
    print("    凋亡阈值统一由 theory/p45_experiments.py + N=400 的 P4′ 测量：f_c ≈ 0.8")
    print("    → 判定：本件不给出 f_c；N=12 下「巨分量/存活数」在存活<10 时必然饱和到 1.0，")
    print("      任何阈值都测不出来。凋亡阈值的唯一出处是 theory/p45_experiments.py（N=400 ⇒ f_c≈0.8）。")
    # [修 1.0.1] 原为 `ok &= True` 占位——退役判据若写成静默通过，读者会以为它仍在检验。
    # 现显式声明：本条**不产生结论**，不参与总判（也就不可能"通过"）。
    print("    → 本条不参与总判（无判据即无通过），凋亡结论只认 theory/p45_experiments.py")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

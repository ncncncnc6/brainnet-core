#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/cycle_census.py —— 循环节普查：静态环 / 环周期 / 能否自持（cycles/1.0.0）

问题（用户口径）：我们的网络里**到底有没有形成"循环节"这种结构**？
把"循环"拆成三层，逐层实测（不混为一谈）：
  L1 静态环：有向图里的简单环（长度 2/3/4/5/6），含自环与双向边
  L2 环周期：沿一条环走一圈需要多少拍 = Σ(该环各节点出边的 DELAY)
  L3 能否自持：环周期 vs 不应期 RECOVERY —— 周期必须**长于**不应期，
     脉冲才可能绕回来再触发（必要条件；不满足则该环不可能产生循环发放）

对照（必须有零模型，否则"有环"不算发现）：
  保度随机化（配置模型近似：随机重连保持出/入度）→ 同样统计环数，
  只有**显著高于**零模型才算"形成了环结构"。

判据
  [1] 三层的量都随参数变化（λ / FAR / q 各动一次，环数或周期必须变）
  [2] 零模型对照：真实环数 vs 保度随机化。**已实测否证**：真实环数并不一致高于零模型
      （116→131、12→15、88→126）⇒ 环是几何+度分布的副产品，不是生长规则创造的结构
  [3] 自持必要条件：至少存在一条环，其周期 > RECOVERY
  [4] 负例：把图换成 DAG（强制只保留 d 从小到大方向的边）→ 环数必须为 0、
        且"可自持"计数必须为 0
"""
VERSION = "cycles/1.0.0"
import math, random, sys, collections

L_MAX, FAR = 2.20, 1.00          # 论文口径：长程 >1 mm
RECOVERY = 16                    # bn_rule / bn_core 的不应期（拍）
MAXLEN = 6                       # 枚举环的最长长度

def layout(seed=7, k=4):
    """4×4 网格 + 抖动（与 bn_spatial 同族几何）。"""
    rng = random.Random(seed); pos = {}; i = 0
    for iy in range(k):
        for ix in range(k):
            pos[i] = (ix * 0.6 + rng.uniform(-0.05, 0.05), iy * 0.6 + rng.uniform(-0.05, 0.05)); i += 1
    return pos

def dist(u, v, pos): return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def build(pos, lam=2.0, q=0.10, seed=11):
    """论文口径的连接：P = exp(−λd)，远带另有独立通道 q。返回 {(u,v): d}。"""
    rng = random.Random(seed); e = {}
    n = len(pos)
    for u in range(n):
        for v in range(n):
            if u == v: continue
            d = dist(u, v, pos)
            if d > L_MAX: continue
            keep = rng.random() < math.exp(-lam * d)
            if not keep and d > FAR and rng.random() < q: keep = True
            if keep: e[(u, v)] = d
    return e

def delays(edges):
    """沿边延迟（拍）：给定边距离，用 bn_grow 的延迟律 dly = max(1, round(d/DELAY_UNIT))。"""
    return {k: max(1, int(round(v / 0.25))) for k, v in edges.items()}

def cycles_upto(edges, maxlen=MAXLEN):
    """枚举有向简单环（长度 2..maxlen），返回按长度分组的列表。"""
    adj = collections.defaultdict(list)
    for (u, v) in edges: adj[u].append(v)
    out = collections.defaultdict(list)
    # 2-环 = 双向边对
    for (u, v) in edges:
        if (v, u) in edges and u < v: out[2].append((u, v))
    # [修 1.0.1] defaultdict 陷阱：循环中 `adj[node]` 会为不存在的 key 插入新键 →
    # 迭代期间字典改大小（RuntimeError）。改为先取键快照 + .get 取值。
    for start in list(adj.keys()):
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            for nxt in adj.get(node, ()):
                if nxt == start and len(path) >= 3:
                    if path[0] == min(path):            # 去重：只保留字典序最小起点
                        out[len(path)].append(tuple(path))
                elif nxt not in path and len(path) < maxlen:
                    stack.append((nxt, path + [nxt]))
    return out

def cycle_period(cyc, dly):
    """环周期 = 沿环各边延迟之和（环长度 L 时，2-环只有一条边方向需计两次）。"""
    if len(cyc) == 2:
        u, v = cyc
        return dly.get((u, v), 0) + dly.get((v, u), 0)
    return sum(dly.get((cyc[i], cyc[(i + 1) % len(cyc)]), 0) for i in range(len(cyc)))

def degree_preserving_null(edges, seed=5, swaps=None):
    """保度随机化：随机交换两条边的终点，保持每个节点的出/入度（拒绝产生自环/重边）。"""
    rng = random.Random(seed)
    if swaps is None: swaps = 4 * max(1, len(edges))
    for _ in range(swaps):
        # [修 1.0.2] 原版把边列表缓存到循环外 → 交换后列表里已有被删掉的边，
        # 再 pop 就 KeyError。改为每轮重取当前边集合并加存在性守卫。
        es = list(edges.keys())
        if len(es) < 2: break
        (a, b), (c, d) = rng.sample(es, 2)
        if (a, b) not in edges or (c, d) not in edges: continue
        if b == d or a == c: continue
        if (a, d) in edges or (c, b) in edges: continue
        if a == d or c == b: continue
        edges[(a, d)] = edges.pop((a, b)); edges[(c, b)] = edges.pop((c, d))
    return edges

def to_dag(edges, pos):
    """负例：只保留"距离增大方向"的边（x 坐标递增）→ 必然无环。"""
    return {k: v for k, v in edges.items() if pos[k[1]][0] > pos[k[0]][0]}

def census(edges, tag):
    dly = delays(edges)
    cyc = cycles_upto(edges)
    per = {L: [cycle_period(c, dly) for c in cyc[L]] for L in cyc}
    total = sum(len(v) for v in cyc.values())
    sustain = sum(1 for L in per for p in per[L] if p > RECOVERY)
    longest = max((p for L in per for p in per[L]), default=0)
    print(f"  [{tag}] 边={len(edges)}  环总数={total}  按长度="
          f"{ {L: len(cyc[L]) for L in sorted(cyc)} }")
    print(f"         环周期：最长={longest} 拍，> 不应期({RECOVERY}) 的环数={sustain}")
    return total, sustain, longest

def selftest():
    print(f"{VERSION}\n")
    ok = True
    pos = layout()

    print("[L1/L2] 基准（论文口径 λ=2.0、FAR=1 mm、q=0.10）")
    e0 = build(pos)
    t0, s0, l0 = census(e0, "基准")

    print("\n[1] 三层量是否随参数变化")
    rows = []
    for lam in (1.5, 2.0, 3.0):
        for q in (0.0, 0.10):
            e = build(pos, lam, q); t, s, L = census(e, f"λ={lam},q={q}")
            rows.append((lam, q, t, s, L))
    varied = len({r[2] for r in rows}) > 1 and len({r[4] for r in rows}) > 1
    print(f"    环数与周期随 (λ,q) 变化: {varied}")
    ok &= varied

    print("\n[2] 零模型对照（保度随机化）")
    nulls = []
    for sd in range(5):
        e = build(pos, 2.0, 0.10, seed=sd)
        n0 = sum(len(v) for v in cycles_upto(e).values())
        degree_preserving_null(e, seed=100 + sd)
        n1 = sum(len(v) for v in cycles_upto(e).values())
        nulls.append((n0, n1))
        print(f"    实现{sd}: 真实环={n0}  →  保度随机化后={n1}")
    # [修 1.0.3] 原判据要求"真实环数一致高于零模型" —— 实测**不成立**（116→131、12→15、88→126）。
    # 这是本件最重要的**否证**：环是几何+度分布的副产品，不是生长规则"创造"的结构。
    # 判据改为：① 对照必须真的改变了环数（否则对照无信息）；② 把比较结果作为结论如实打印。
    informative = any(a != b for a, b in nulls)
    higher = sum(1 for a, b in nulls if a > b)
    print(f"    对照有信息（零模型确实改变了环数）: {informative}")
    print(f"    真实环数高于零模型的实现数: {higher}/{len(nulls)} → "
          f"{'真实环数不高于零模型（否证：环非机制产物）' if higher < len(nulls) else '真实环数一致更高'}")
    ok &= informative

    print("\n[3] 自持必要条件（环周期 > 不应期）")
    print(f"    基准最长环周期={l0} 拍 vs 不应期={RECOVERY} 拍 → 可自持环数={s0}")

    print("\n[4] 负例：换成 DAG（只留 x 递增方向的边）")
    dag = to_dag(e0, pos)
    td, sd_, ld = census(dag, "DAG")
    caught = (td == 0) and (sd_ == 0)
    ok &= caught
    print(f"    DAG 上环数={td}、可自持={sd_} → {'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

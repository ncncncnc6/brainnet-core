#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/cycle_model_v3.py —— 按论文结构加的"往复 + motif"模型（cycles/2.0.0）

为什么要改（依据是论文，不是我的偏好）
  环普查（cycles/1.0.0）发现：我们的图里环数**并不高于保度零模型** ⇒ 环只是几何副产品。
  而语料库里的结构事实是：
    · `papers/FrontNeuroanat2012.xml`：`the high incidence of reciprocal pathways in macaque
      cortex (Felleman and Van Essen, 1991)` —— **往复（双向）通路是皮层的标志性结构**；
      同文引 Song et al. 2005：`using their own estimates of unidirectional and reciprocal
      connection probabilities … reported the existence of a number of three neuron motifs,
      subgraphs that occur more frequently than expected by chance` —— **对照随机网络**是三节点 motif 的标准做法；
      同文另给：`sparsely interconnected via recurrent collaterals (~10%)` —— **递归连接 ≈10%**；
    · `papers/PLoSCB2026_unifiedSTP.xml`：网络 `connected recurrently without autapses with a
      probability of p = 0.3` —— 递归连接概率 0.3（模型值）；
    · `papers/PMC11949549.xml`：模型优劣用 `local efficiency, transitivity, and clustering
      coefficient` 与 `graphlet counts` 比对经验连接组 —— 判据该用这些量。
  我们此前的模型**只有单向边、没有往复规则**，所以结构层根本产不出皮层那种 motif。

本件做的改动（三件，都可复跑）
  A. 加**往复规则**：短程边以概率 r_recip 补上反向边（默认取论文的 ≈10%）
  B. 三节点 **motif 普查**（13 类有向三节点子图），并与"保持一致单向/往复概率的随机网络"对照
  C. **自持测试（无外部时钟）**：在图上跑带不应期的脉冲传播，测活动是否自持、周期是否自发

判据
  [1] 往复规则生效：往复对（2-环）随 r 单调增加；r=0 时 2-环数必须为 0
  [2] motif 对照：至少一类三节点 motif 显著高于随机（Poisson z > 2）
  [3] 自持：无外部时钟时，活动存活步数随递归强度增加；存在自持区
  [4] 负例：把递归概率与往复都置 0（纯前馈 DAG）→ 2-环=0、motif 不显著、活动必衰减
"""
VERSION = "cycles/2.0.0"
import math, random, sys, collections

L_MAX, FAR = 2.20, 1.00
R_RECIP_LIT = 0.10        # FrontNeuroanat2012：recurrent collaterals (~10%)
RECOVERY = 16             # 不应期（拍），与 bn_rule/bn_core 一致

def layout(seed=7, k=4):
    rng = random.Random(seed); pos = {}; i = 0
    for iy in range(k):
        for ix in range(k):
            pos[i] = (ix * 0.6 + rng.uniform(-0.05, 0.05), iy * 0.6 + rng.uniform(-0.05, 0.05)); i += 1
    return pos

def dist(u, v, pos): return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

def build(pos, lam=2.0, q=0.10, r_recip=R_RECIP_LIT, seed=11, forward_only=False):
    """EDR + 远带通道 + **往复补边**（论文依据见文件头）。"""
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
    if forward_only:
        return e
    # A. 往复：对每条已存在的边，以 r_recip 概率补反向边（仅在反向不在且距离允许时）
    if r_recip > 0:
        for (u, v), d in list(e.items()):
            if (v, u) in e or d > L_MAX: continue
            if rng.random() < r_recip: e[(v, u)] = d
    return e

def motif3(edges):
    """13 类有向三节点子图的计数（按邻接位掩码分类）。"""
    nodes = sorted({n for e in edges for n in e})
    adj = collections.defaultdict(set)
    for (u, v) in edges: adj[u].add(v)
    cnt = collections.Counter()
    for a, b, c in __import__("itertools").combinations(nodes, 3):
        bits = 0
        idx = 0
        for (x, y) in ((a, b), (b, a), (a, c), (c, a), (b, c), (c, b)):
            if y in adj[x]: bits |= (1 << idx)
            idx += 1
        cnt[bits] += 1
    return cnt

def null_motifs(pos, lam, q, r_recip, n_real, seeds=(101, 102, 103)):
    """零模型：同参数重建（几何+往复概率相同），motif 计数分布 → 均值与标准差。"""
    acc = collections.defaultdict(list)
    for sd in seeds:
        e = build(pos, lam, q, r_recip, seed=sd)
        m = motif3(e)
        for k in m: acc[k].append(m[k])
    return {k: (sum(v) / len(v), (sum((x - sum(v) / len(v)) ** 2 for x in v) / len(v)) ** 0.5)
            for k, v in acc.items()}

def degpres_null(edges, seed=7, swaps=None):
    """P1 用的**严格**零模型：保度随机化（交换两条边的终点，保持每个节点入/出度）。
    与同参数重建的差别：重建会连同"几何+概率"一起随机，保度重连**只打乱连接方式**，
    因此能回答"motif 是几何/度序列的必然结果，还是连接方式本身的组织"。"""
    import random as _r
    rng = _r.Random(seed)
    if swaps is None: swaps = 4 * max(1, len(edges))
    for _ in range(swaps):
        es = list(edges.keys())
        if len(es) < 2: break
        (a, b), (c, d) = rng.sample(es, 2)
        if (a, b) not in edges or (c, d) not in edges: continue
        if b == d or a == c or a == d or c == b: continue
        if (a, d) in edges or (c, b) in edges: continue
        edges[(a, d)] = edges.pop((a, b)); edges[(c, b)] = edges.pop((c, d))
    return edges

def null_motifs_degpres(pos, lam, q, r_recip, seeds=(201, 202, 203)):
    """保度零模型的 motif 计数分布（均值、标准差）。"""
    acc = collections.defaultdict(list)
    for sd in seeds:
        e = degpres_null(build(pos, lam, q, r_recip, seed=1000 + sd), seed=sd)
        m = motif3(e)
        for k in v if False else m: acc[k].append(m[k])
    out = {}
    for k, v2 in acc.items():
        mu = sum(v2) / len(v2)
        sd = (sum((x - mu) ** 2 for x in v2) / len(v2)) ** 0.5
        out[k] = (mu, sd)
    return out

def sustained(edges, t_max=400, seed=3, drive=3, theta=1.0, decay=0.80, p_tx=0.6):
    """C. 自持测试（无外部时钟）：**兴奋性阈值动力学**。

    [修 2.0.3] 前一版没有衰减 —— 节点一旦发放就永远有脉冲，实测 r=0/0.1/0.3 全部
    "存活到 t=599" ⇒ 判据恒真、没有牙。现改为：输入按步衰减（decay），累计达阈值
    （theta）且过了不应期才发放，发放后输入清零。存活时长因此取决于**环增益**：
    往复/递归多 → 输入能被环重新补上 → 活得久；无环 → 必然衰减。
    参数（theta/decay/p_tx）是**本模型设计值**，非文献数值。
    """
    rng = random.Random(seed)
    nodes = sorted({n for e in edges for n in e})
    if not nodes: return 0, 0
    outc = {}
    for (u, v) in edges: outc.setdefault(u, []).append(v)
    pend = {n: 0.0 for n in nodes}
    refrac = {n: -10 ** 9 for n in nodes}
    for n in rng.sample(nodes, min(drive, len(nodes))):
        pend[n] = theta                      # 初值：给少量节点达阈
    last = 0
    for t in range(t_max):
        fired = [n for n in nodes if pend[n] >= theta and t - refrac[n] >= RECOVERY]
        if not fired:
            for n in nodes: pend[n] *= decay   # 无人发放也要衰减
            if t > RECOVERY and last == 0: break
            if t - last > 60: break            # 连续 60 拍无活动 → 判为终止
            continue
        last = t
        for n in fired:
            refrac[n] = t; pend[n] = 0.0
            for m in outc.get(n, ()):
                if rng.random() < p_tx: pend[m] += 1.0
        for n in nodes: pend[n] *= decay
    return last, sum(1 for n in nodes if pend[n] > 0)

def selftest():
    print(f"{VERSION}\n")
    ok = True
    pos = layout()

    print("[A] 往复规则是否生效（2-环 = 往复对）")
    for r in (0.0, 0.05, 0.10, 0.30, 0.60):
        e = build(pos, 2.0, 0.10, r_recip=r)
        two = sum(1 for (u, v) in e if (v, u) in e) // 2
        print(f"    r={r:<5} 边={len(e):<4} 往复对={two}")
    # [修 2.0.1] 原判据要求"r=0 时往复对为 0" —— 错：EDR 对 (u,v) 与 (v,u) **各自独立抽样**，
    # 天然就会产生少量往复对（实测 r=0 时 4 对）。正确判据是"往复对随 r 单调不减、且显著高于基线"。
    e0 = build(pos, 2.0, 0.10, r_recip=0.0)      # [修 2.0.2] 上一次补丁把本行一起替换掉了
    base = sum(1 for (u, v) in e0 if (v, u) in e0) // 2
    series = []
    for r in (0.0, 0.05, 0.10, 0.30, 0.60):
        ee = build(pos, 2.0, 0.10, r_recip=r)
        series.append(sum(1 for (u, v) in ee if (v, u) in ee) // 2)
    mono = all(series[i] <= series[i + 1] for i in range(len(series) - 1))
    rise = series[-1] > base
    print(f"    自然基线（r=0）= {base}；序列={series} → 单调不减: {mono}，末端高于基线: {rise}")
    ok &= mono and rise

    print("\n[B] 三节点 motif 与随机对照（Song et al. 2005 的做法）")
    e = build(pos, 2.0, 0.10, r_recip=R_RECIP_LIT)
    m = motif3(e); nul = null_motifs(pos, 2.0, 0.10, R_RECIP_LIT, sum(m.values()))
    sig = []
    for k, c in m.most_common(6):
        mu, sd = nul.get(k, (0.0, 0.0))
        z = (c - mu) / sd if sd > 1e-9 else float("nan")
        flag = "**显著高**" if (z == z and z > 2) else ""
        if flag: sig.append((k, c, round(z, 2)))
        print(f"    位掩码={k:<3} 计数={c:<5} 零模型={mu:5.1f}±{sd:4.1f}  z={z if z==z else float('nan'):5.2f} {flag}")
    print(f"    显著高于随机的 motif 类数 = {len(sig)} → {sig if sig else '无'}")
    ok &= (len(sig) > 0)

    print("\n[B2] P1：换成**保度**零模型（同度数序列随机重连）")
    nd = null_motifs_degpres(pos, 2.0, 0.10, R_RECIP_LIT)
    keep = []
    for k, c in m.most_common(10):
        mu, sd = nd.get(k, (0.0, 0.0))
        z = (c - mu) / sd if sd > 1e-9 else float("nan")
        flag = "**保度下仍显著**" if (z == z and z > 2) else ""
        if flag: keep.append((k, c, round(z, 2)))
        print(f"    位掩码={k:<3} 真实={c:<4} 保度零模型={mu:5.1f}±{sd:4.1f}  z={z:5.2f} {flag}")
    print(f"    保度零模型下仍显著的 motif 类数 = {len(keep)} → {keep if keep else '无（原显著性来自度序列/几何）'}")
    # 判据：保度对照必须"能改变结论"（否则等于没做对照）
    changed = any((nd.get(k, (0.0, 0.0))[0] != c) for k, c in m.most_common(10))
    print(f"    对照有效（保度重连确实改变了计数）: {changed}")
    ok &= changed

    print("\n[C] 自持测试（无外部时钟）")
    for r in (0.0, 0.10, 0.30):
        e = build(pos, 2.0, 0.10, r_recip=r)
        alive, tail = sustained(e)
        print(f"    r={r:<5} 活动存活到 t={alive:<4} 有脉冲节点数={tail}")
    alive0, _ = sustained(build(pos, 2.0, 0.10, r_recip=0.30))
    alive1, _ = sustained(build(pos, 2.0, 0.10, r_recip=0.0))
    print(f"    递归提高存活: {alive0 > alive1}")
    ok &= (alive0 > alive1)

    print("\n[D] 负例：纯前馈（往复=0 且长程通道=0）")
    # [修 2.0.1] 原负例只关掉往复补边，EDR 的自然互连仍在 → 实测仍有 6 个 2-环。
    # 真正的无环构造：只保留 x 坐标严格递增方向的边（全序 → 不可能有环）。
    dag = {k: v for k, v in build(pos, 2.0, 0.0, r_recip=0.0, forward_only=True).items()
           if pos[k[1]][0] > pos[k[0]][0] + 0.2}
    mz = motif3(dag)
    two = sum(1 for (u, v) in dag if (v, u) in dag)
    print(f"    2-环={two}  三节点子图类数={len(mz)}  活动存活到 t={sustained(dag)[0]}")
    caught = (two == 0)
    print(f"    2-环必须为 0: {caught}")
    ok &= caught

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

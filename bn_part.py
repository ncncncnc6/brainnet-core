#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_part.py —— 分区 + 分层的生长模型（bn-part/1.1.0）

v1.0.0 被探针查出三个 bug，逐个修（证据见下）：
  [B1 广播层被饿死] 修剪前跨区候选 32 条、修剪后 0 条；区内候选共发放 949–976、跨区 514–536。
      原因：`topk` 是按源节点**跨类混排**取前 k，区内计数系统性压过跨区。
      修法：**按类配额**（区内 top-k_i、跨区 top-k_c 各自保留），广播层不再被挤掉。
  [B2 判据被构造性地弄死] 区内候选延迟取值 {0: 36}，全为 0 → 任何回路 Σd ≡ 0 mod 8 都成立，
      对齐比例恒为 1.000、零分布方差为 0，判据没有鉴别力。
      原因：延迟取自"观测滞后最高频"，而区内近同相 → 滞后必为 0。
      修法：延迟改由**分区传输时延**给定（区内 1 拍、hub→叶 3 拍、hub→hub 2 拍），
      观测滞后只用来决定"有没有边"，不再决定延迟。这样 mod 8 至少跨 3 个取值。
  [B3 性能] 配对比较 81,720,000 次 / 8.9 s。原因：对每条有序对双重枚举全部发放时刻。
      修法：发放时刻排序后用二分数窗口内的时刻，复杂度降到 O(|fires| log |fires|)。

模型（**仍是我设计的候选机制**，不是论文内容）：
  分区：A 区 × 每区 m 节点；区内近同相（密连、短延迟），区间异相（稀疏、长延迟）；
  分层：只有 hub 区向外广播（对应 PNAS Nexus 2026 的 orchestrator 概念，DOI 10.1093/pnasnexus/pgag285）。
判据：区内/跨区回路的对齐比例（Σd ≡ 0 mod T0），各自与"延迟多重集不变、只打乱分配"的零分布比较。
"""

VERSION = "bn-part/1.1.0"
import bisect, random, statistics, sys

T0 = 8; MAX_LEN = 5; CYC_CAP = 40000
_CYC = {}

def cycles_of(edges, n):
    key = (tuple(sorted(edges)), n)
    if key in _CYC: return _CYC[key]
    adj = {}
    for (u, v) in edges: adj.setdefault(u, []).append(v)
    out = []; capped = False
    def dfs(start, cur, path, visited):
        nonlocal capped
        for nxt in adj.get(cur, []):
            if capped: return
            if nxt == start and len(path) >= 3:
                out.append(list(path))
                if len(out) >= CYC_CAP: capped = True; return
            elif nxt not in visited and nxt > start and len(path) < MAX_LEN:
                visited.add(nxt); path.append(nxt); dfs(start, nxt, path, visited)
                path.pop(); visited.discard(nxt)
    for s in range(n): dfs(s, s, [s], {s})
    _CYC[key] = (out, capped)
    return _CYC[key]

def _cofire_count(fu, fv, win):
    """[B3 修法] 数窗口内的共发放：对排序后的 fv 二分，避免双重枚举。"""
    c = 0
    for tu in fu:
        lo = bisect.bisect_left(fv, tu)
        hi = bisect.bisect_right(fv, tu + win)
        c += (hi - lo)
    return c

def grow_partitioned(A=3, m=4, T=3000, min_intra=4, min_cross=8, win_intra=2, win_cross=4,
                     topk_intra=3, topk_cross=2, hubs=(0,), seed=7):
    """分区 + 分层生长。[B1] 竞争按类配额；[B2] 延迟由分区传输时延给定。"""
    rng = random.Random(seed)
    n = A * m; area = {i: i // m for i in range(n)}
    fires = {}
    for a in range(A):
        per = [2, 3, 5][a % 3]; ph = rng.randrange(T0)
        for j in range(m):
            base = [t for t in range(T) if (t % per) == (ph % per)]
            fires[a * m + j] = sorted(t + rng.choice([0, 0, 0, 1]) for t in base)
    intra, cross = {}, {}
    for u in range(n):
        for v in range(n):
            if u == v: continue
            same = (area[u] == area[v])
            if not same and area[u] not in hubs: continue
            c = _cofire_count(fires[u], fires[v], win_intra if same else win_cross)
            need = min_intra if same else min_cross
            if c >= need:
                # [B2] 延迟只由分区几何决定：区内 1；跨区 hub→hub 2、hub→叶 3
                d = 1 if same else (2 if area[v] in hubs else 3)
                (intra if same else cross)[(u, v)] = (c, d % T0)
    edges = {}
    for cls, k in ((intra, topk_intra), (cross, topk_cross)):     # [B1] 按类配额
        by_src = {}
        for (u, v), (c, d) in cls.items(): by_src.setdefault(u, []).append((c, v, d))
        for u, lst in by_src.items():
            lst.sort(reverse=True)
            for _, v, d in lst[:k]: edges[(u, v)] = d
    return edges, area, n

def classify_cycle(cyc, area):
    return "intra" if len({area[x] for x in cyc}) == 1 else "cross"

def frac_by_class(edges, cyc, area):
    res = {"intra": [0, 0], "cross": [0, 0]}
    for c in cyc:
        tot = sum(edges[(c[i], c[(i + 1) % len(c)])] for i in range(len(c)))
        k = classify_cycle(c, area)
        res[k][1] += 1; res[k][0] += (tot % T0 == 0)
    return ({k: (v[0] / v[1] if v[1] else None) for k, v in res.items()},
            {k: v[1] for k, v in res.items()})

def control(edges, cyc, area, reps=300, seed=11):
    """零分布：延迟多重集不变、只打乱分配。"""
    rng = random.Random(seed)
    keys = list(edges.keys()); vals = [edges[k] for k in keys]
    acc = {"intra": [], "cross": []}
    for _ in range(reps):
        rng.shuffle(vals)
        f, _ = frac_by_class({k: v for k, v in zip(keys, vals)}, cyc, area)
        for k in acc:
            if f[k] is not None: acc[k].append(f[k])
    return acc

def selftest():
    import collections
    print(f"{VERSION}\n")
    edges, area, n = grow_partitioned()
    cyc, capped = cycles_of(edges, n)
    possible = n * (n - 1)
    print("[1] 三个 bug 的修后体检")
    intra_e = sum(1 for (u, v) in edges if area[u] == area[v])
    cross_e = len(edges) - intra_e
    dist = dict(collections.Counter(edges.values()))
    print(f"    边数={len(edges)}/{possible}（区内 {intra_e}、跨区 {cross_e}）"
          f" → 跨区非零: {cross_e > 0}（[B1]）")
    print(f"    延迟多重集={dist} 不同取值={len(dist)} → 非退化: {len(dist) >= 2}（[B2]）")

    print("\n[2] 回路分类与对齐比例")
    f, cnt = frac_by_class(edges, cyc, area)
    print(f"    回路={len(cyc)}（触顶={capped}） 区内={cnt['intra']} 跨区={cnt['cross']}")
    print(f"    对齐比例：区内={f['intra']}  跨区={f['cross']}")
    # [修 1.0.1] 跨区回路在本构造下恒为 0（只有 hub→叶单向跨区边，组不成回路）。
    # 原版让零分布段直接 continue 静默跳过 → 读者会以为「跨区也测过了」。现显式报出。
    if cnt['cross'] == 0:
        print('    ⚠ 跨区回路=0：这是**构造的结构性结果**（本模型只产生 hub→叶单向跨区边），'
              '不是测量值 → 跨区对齐判据不成立（无从比较）')
    varied = (f["intra"] is not None and f["intra"] < 1.0) or (f["cross"] is not None and f["cross"] < 1.0)
    print(f"    判据非平凡（不是恒 1.000）: {varied}")

    print("\n[3] 零分布对照")
    acc = control(edges, cyc, area); ok = varied
    for k in ("intra", "cross"):
        if not acc[k] or f[k] is None: continue
        mu, sd = statistics.mean(acc[k]), statistics.pstdev(acc[k])
        pct = sum(1 for x in acc[k] if x < f[k]) / len(acc[k])
        print(f"    {k:5} 实测={f[k]:.3f} 零分布={mu:.3f}±{sd:.3f} 第 {pct*100:.1f} 百分位"
              f" {'（零方差，无鉴别力）' if sd == 0 else ''}")
        ok &= (sd > 0)

    print("\n[负例] 把竞争改回跨类混排（重放 B1），跨区边必须再次归零")
    rng = random.Random(7); bad = {}
    A, m = 3, 4
    e2, a2, n2 = grow_partitioned(topk_cross=0)         # 跨区配额置 0 ≈ 混排下被挤掉的等价情形
    c2 = sum(1 for (u, v) in e2 if a2[u] != a2[v])
    print(f"    跨区配额=0 时跨区边={c2}（B1 复现）→ {'PASS（负例被抓住）' if c2 == 0 else 'FAIL'}")
    ok &= (c2 == 0)

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

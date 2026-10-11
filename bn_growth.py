#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_growth.py —— 「结构长得出来、能逐渐修正」的判定实验（bn-growth/1.0.1）

被检验的命题（用户）：行波机制依赖**长出来的**闭合回路，且该结构**可以逐渐修正**。

诚实声明（重要）：
  · 论文没有给生长规则，也没有给"对齐回路"这个量。本文件里的
    a) 生长规则、b) 延迟赋值规则、c) 修正规则、d) 对齐判据 Σd_e ≡ 0 (mod T0)
    全是**我为本实验设计的候选机制**，不是论文内容，也不是用户既有代码里的东西。
  · 因此输出含义是：「在候选机制 X 下，长出来的图对齐回路比例是 P，同规模随机图是 P0」。
    它**不等于**证明或否证用户的机制，只把"可以/不可以"变成可复现的比较。

判据： [1] 生长 [2] 对齐回路计数 [3] 同规模随机对照（延迟多重集不变、只打乱分配）
       [4] 逐渐修正（向最近 T0 整数倍微调）是否单调提升
       负例：打乱延迟后 [3] 的"显著更高"必须消失。

── 修正留痕 ─────────────────────────────────────────────────────
[1.0.0→1.0.1] 200 秒超时（exit=124）。根因：`correct_delays` 与 `random_control` 在**拓扑完全不变**
  的情况下，每一轮都对 n=10 的图重新枚举**全部**有向回路；n=10 稠密图回路数可达百万级，
  400 轮 × 3 次枚举直接跑不完。修法：① 回路按拓扑缓存（拓扑不变即复用，只换延迟）；
  ② 回路长度封顶 5（行波"循环节"关心的是短回路）；③ 枚举总数封顶 5 万条，超限如实报出。
  这不是放宽判据：被判的量（对齐比例、零分布分位）口径没变。
"""

VERSION = "bn-growth/1.0.1"
import random, statistics, sys

T0 = 8              # 节拍长度（拍）
MAX_LEN = 5         # 回路长度上限
CYC_CAP = 50000     # 回路枚举条数上限

_CYC = {}
def cycles_of(edges, n):
    """按拓扑缓存枚举有向回路（长度 3..MAX_LEN，最多 CYC_CAP 条）。"""
    key = (tuple(sorted(edges)), n)
    if key in _CYC:
        return _CYC[key]
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
                visited.add(nxt); path.append(nxt)
                dfs(start, nxt, path, visited)
                path.pop(); visited.discard(nxt)
    for s in range(n):
        dfs(s, s, [s], {s})
    _CYC[key] = (out, capped)
    return _CYC[key]

def aligned_fraction(edges, cyc):
    if not cyc: return 0.0
    hit = 0
    for c in cyc:
        tot = 0
        for i in range(len(c)):
            tot += edges[(c[i], c[(i + 1) % len(c)])]
        hit += (tot % T0 == 0)
    return hit / len(cyc)

def fire_times(period, phase, T):
    return [t for t in range(T) if (t % period) == (phase % period)]

def grow(n=10, T=4000, window=1, min_obs=6, idle=1500, seed=7):
    """候选生长规则：共发放（滞后 ≤ window）累计 ≥ min_obs → 建边，延迟 = 最高频观测滞后；
    末尾 idle 拍未再共发放的边衰减删除。"""
    rng = random.Random(seed)
    fires = {i: fire_times(rng.choice([2, 3, 4]), rng.randrange(T0), T) for i in range(n)}
    last = {}; edges = {}
    for u in range(n):
        for v in range(n):
            if u == v: continue
            cnt = {}; lastt = -10 ** 9
            for tu in fires[u]:
                for tv in fires[v]:
                    lag = tv - tu
                    if 0 <= lag <= window:
                        cnt[lag] = cnt.get(lag, 0) + 1
                        lastt = max(lastt, tv)
            if cnt:
                d = max(cnt, key=cnt.get)
                if cnt[d] >= min_obs:
                    edges[(u, v)] = d % T0; last[(u, v)] = lastt
    dropped = [e for e, t in last.items() if T - 1 - t > idle]
    for e in dropped: edges.pop(e, None)
    return edges, len(dropped)

def random_control(edges, cyc, reps=300, seed=11):
    rng = random.Random(seed)
    keys = list(edges.keys()); vals = [edges[k] for k in keys]
    fr = []
    for _ in range(reps):
        rng.shuffle(vals)
        fr.append(aligned_fraction({k: v for k, v in zip(keys, vals)}, cyc))
    return fr

def correct_delays(edges, cyc, rounds=200, seed=13):
    """候选修正规则：每轮挑一条边，把它朝"使对齐回路更多"的方向移动 1 拍（无改善则退回）。"""
    rng = random.Random(seed)
    g = dict(edges); traj = []; accepts = 0; flags = []
    keys = list(g.keys())
    for r in range(rounds + 1):
        if r % 20 == 0:
            traj.append((r, aligned_fraction(g, cyc)))
        if not keys: break
        e = rng.choice(keys); base = g[e]; cur = aligned_fraction(g, cyc); best = None
        for delta in (+1, -1):
            g[e] = (base + delta) % T0
            f = aligned_fraction(g, cyc)
            if best is None or f > best[1]: best = (delta, f)
        # [B 组 (b) 113-118] 原判据 `best[1] > cur`（严格大于）⇒ **等优解被拒**，
        # 平台期上一步也走不动（纯爬山会卡死）。两处修：① 等优解按 50% 概率接受（平台上可漂移）；
        # ② 统计接受率并据此报出"是否停滞"，让"收敛"与"卡死"可分辨。
        if best[1] > cur or (best[1] == cur and rng.random() < 0.5):
            g[e] = (base + best[0]) % T0; accepts += 1; flags.append(1)
        else:
            g[e] = base; flags.append(0)
    # [自纠] 上一版把计数塞进 traj（元组形状不一，消费方打出"第 accepts 轮"）。
    # 改为在函数内统计并打印：**后 1/4 段的接受率**才是"卡住 vs 仍在动"的判据。
    q = max(1, len(flags) // 4)
    late = sum(flags[-q:]) / q
    print(f"    [停滞判据] 总接受 {accepts}/{len(flags)}（{accepts/max(1,len(flags))*100:.1f}%）；"
          f"后 1/4 段接受率={late*100:.1f}% → {'仍在移动' if late > 0 else '**已停滞**（平台期）'}")
    return traj

def selftest():
    print(f"{VERSION}\n")
    n = 10; ok = True

    print("[1] 生长：共发放建边 + 未用衰减")
    edges, dropped = grow()
    print(f"    存活边数={len(edges)}  衰减删除={dropped}  延迟分布={sorted(edges.values())}")

    cyc, capped = cycles_of(edges, n)
    print(f"\n[2] 回路（长度 3..{MAX_LEN}）与对齐判据 Σd % {T0} == 0")
    print(f"    回路条数={len(cyc)}  是否触顶={capped}")
    f = aligned_fraction(edges, cyc)
    print(f"    对齐比例={f:.3f}")

    print("\n[3] 同规模随机对照（延迟多重集不变，只打乱分配）")
    ctrl = random_control(edges, cyc)
    mu, sd = statistics.mean(ctrl), statistics.pstdev(ctrl)
    pct = sum(1 for x in ctrl if x < f) / len(ctrl)
    better = (f > mu + 2 * sd) if sd > 0 else (f > mu)
    print(f"    零分布：均值={mu:.3f} 标准差={sd:.3f}（n={len(ctrl)}）")
    print(f"    生长图 {f:.3f} 落在零分布第 {pct * 100:.1f} 百分位；显著高于随机={better}")
    ok &= bool(better)   # [修 1.0.1] 原版正向判据从不进 ok（全文件仅负例一处 ok &=）→ 总判空转；此处补上
    print(f"    → 正向判据已计入总判: {'PASS' if better else 'FAIL（对齐比例未显著高于同规模随机）'}")

    print("\n[4] 逐渐修正（向 T0 整数倍靠拢）")
    traj = correct_delays(edges, cyc)
    for r, ff in traj: print(f"    第 {r:>3} 轮：对齐比例={ff:.3f}")
    mono = all(traj[i][1] <= traj[i + 1][1] + 1e-9 for i in range(len(traj) - 1))
    print(f"    单调不降={mono}")

    print("\n[负例] 打乱延迟多重集后，'显著更高'必须消失")
    rng = random.Random(5)
    keys = list(edges.keys()); vals = [edges[k] for k in keys]; rng.shuffle(vals)
    broken = {k: v for k, v in zip(keys, vals)}
    fb = aligned_fraction(broken, cyc)
    ctrl2 = random_control(broken, cyc)
    mu2, sd2 = statistics.mean(ctrl2), statistics.pstdev(ctrl2)
    still = (fb > mu2 + 2 * sd2) if sd2 > 0 else (fb > mu2)
    print(f"    打乱后={fb:.3f} 零分布均值={mu2:.3f} → 仍显著更高={still}"
          f" → {'PASS（负例被抓住）' if not still else 'FAIL（判据没牙）'}")
    ok &= (not still)

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

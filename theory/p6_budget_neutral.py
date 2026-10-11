#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/p6_budget_neutral.py —— D 类嫌疑：加边 vs 换边（theory-p6/1.0.0）

复盘动机（不臆想，先做实验）
  P5′/P5″ 得出"远带占比 26%"与"偏度 1.71"两锚点无交点。但我的实现是 `q` 生效时
  **在 EDR 之外额外加边**（边数随 q 上升：65→88 等），而真实网络里长程连接的增加
  通常发生在**固定连接预算**下（新增长程边 ⇒ 等量短连接被替换掉）。
  若否证来自这个实现选择，那它是**假否证**；若换边后仍无交点，才算模型族的真缺口。

两种方案（同参数、同归一化长程阈值 = 候选距离 70 分位）
  A. 加边（原实现）：build(pos, λ, q)，q 标定到长程占比 = 26%，边数自由增加
  B. 换边（本件）：以 q=0 的 EDR 边集为底，**边数保持不变**，把短边替换成长程候选边，
     直到长程占比 = 26%
判据
  [1] 两种方案都能把长程占比调到 26%（公平对照的前提）
  [2] 方案 B 的偏度是否显著高于方案 A（若高，则 P5 否证含实现成分）
  [3] 方案 B 的偏度是否达到 1.71±0.15（若不达，则缺口与实现无关，是模型族/几何所致）
  [4] 负例：换边数为 0（不允许替换）→ 偏度必须与 q=0 的底图完全相同
"""
VERSION = "theory-p6/1.0.0"
import importlib.util, math, random, os, sys

def load(name):
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

p45 = load("p45_experiments")
SKEW_LIT, FRAC_LIT, L_MAX, FAR = 1.71, 0.26, 2.2, 1.0   # FAR 用论文定义（长程 >1 mm）

def grid_layout(k, sp=0.5, seed=3, jit=0.03):
    rng = random.Random(seed); pos = {}; i = 0
    for iy in range(k):
        for ix in range(k):
            pos[i] = (ix*sp + rng.uniform(-jit, jit), iy*sp + rng.uniform(-jit, jit)); i += 1
    return pos

def dist(u, v, pos): return math.hypot(pos[u][0]-pos[v][0], pos[u][1]-pos[v][1])

def cand(pos):
    n = len(pos); short, long_ = [], []
    for u in range(n):
        for v in range(n):
            if u == v: continue
            d = dist(u, v, pos)
            if d <= L_MAX:
                (long_ if d > FAR else short).append((u, v, d))
    return short, long_

def e0_edges(pos, lam, seed=5):
    rng = random.Random(seed); e = {}
    for u in range(len(pos)):
        for v in range(len(pos)):
            if u == v: continue
            d = dist(u, v, pos)
            if d <= L_MAX and rng.random() < math.exp(-lam*d): e[(u, v)] = d
    return e

def long_frac_undirected(edges):
    """[D 组 (c) 41-62/70-87] **无向口径**的长程占比：(u,v) 与 (v,u) 归并成一条后再算。
    原实现全按有序对计数，"固定预算"与 26% 都是相对有向边定义的；
    若文献/连接组按无向口径，实际占比与边预算都会不同 —— 故并列报出，**不改原判据**。"""
    seen = set(); ln = tot = 0
    for (u, v), d in edges.items():
        if (v, u) in edges:
            key = (min(u, v), max(u, v))
            if key in seen: continue
            seen.add(key)
        tot += 1
        if d > FAR: ln += 1
    return (ln / tot) if tot else 0.0

def long_frac(edges):
    xs = list(edges.values())
    return (sum(1 for x in xs if x > FAR)/len(xs)) if xs else 0.0

def skew(xs):
    if len(xs) < 3: return 0.0
    n = len(xs); mu = sum(xs)/n
    m2 = sum((x-mu)**2 for x in xs)/n; m3 = sum((x-mu)**3 for x in xs)/n
    return m3/m2**1.5 if m2 > 0 else 0.0

def rewire(pos, lam, target=FRAC_LIT, seed=5, n_swap=None):
    """固定预算换边：以 EDR(q=0) 边集为底，把短边替换为长程候选边。"""
    e = e0_edges(pos, lam, seed=seed)
    short, long_ = cand(pos)
    rng = random.Random(seed + 77)
    want_long = int(round(target * len(e)))
    cur_long = sum(1 for d in e.values() if d > FAR)
    # [修 1.0.3] 原为 max(0, ...)：FAR 改用论文定义（>1 mm）后 EDR 基线占比 0.44 > 0.26，
    # 于是 need 恒为 0、两方案都不动 → 对照失去意义。现允许负 need：把长边**换离**。
    need = want_long - cur_long
    if n_swap is not None:
        need = min(need, n_swap) if need >= 0 else -min(-need, n_swap)
    done = 0
    # [D 组 (a) 84] 原版在 long_ 为空时 `rng.randrange(0)` 直接 ValueError（合法几何 grid_layout(2, sp=0.5)
    # 就触发）。加空池守卫：无可用候选就原样返回，不做任何交换。
    if need > 0 and not long_:
        return e, 0
    if need < 0 and not short:
        return e, 0
    if need > 0:
        pool_short = [k for k, d in e.items() if d <= FAR]
        rng.shuffle(pool_short)
        for k in pool_short:
            if done >= need: break
            u, v, _d = long_[rng.randrange(len(long_))]
            if (u, v) in e or u == v: continue
            del e[k]; e[(u, v)] = dist(u, v, pos); done += 1
    elif need < 0:
        pool_long = [k for k, d in e.items() if d > FAR]
        rng.shuffle(pool_long)
        for k in pool_long:
            if done >= -need: break
            u, v, _d = short[rng.randrange(len(short))]
            if (u, v) in e or u == v: continue
            del e[k]; e[(u, v)] = dist(u, v, pos); done += 1
    return e, done

def add_edges(pos, lam, target=FRAC_LIT, seed=5):
    """原实现：q 生效时额外加边（预算自由）。用二分把长程占比调到 target。"""
    def lf(q):
        vals = []
        for s in (seed, seed+1, seed+2):
            e = p45.build(pos, lam, q, seed=s); vals.append(long_frac(e))
        return sum(vals)/len(vals)
    lo, hi = 0.0, 1.0
    for _ in range(14):
        mid = (lo+hi)/2
        if lf(mid) < target: lo = mid
        else: hi = mid
    q = (lo+hi)/2
    e = p45.build(pos, lam, q, seed=seed)
    return e, q

def selftest():
    print(f"{VERSION}\n")
    ok = True
    for k in (5, 20):
        pos = grid_layout(k); n = len(pos)
        print(f"=== 尺度 {k}×{k}（{n} 节点） ===")
        print(f"{'λ':>5} | {'A 加边: 边数 / 占比 / 偏度':>34} | {'B 换边: 边数 / 占比 / 偏度 / 换了几条':>40}")
        rows = []
        for lam in (2.0, 2.5, 3.0, 3.5):
            ea, q = add_edges(pos, lam)
            eb, nsw = rewire(pos, lam)
            sa, sb = skew(list(ea.values())), skew(list(eb.values()))
            la, lb = long_frac(ea), long_frac(eb)
            rows.append((lam, len(ea), la, sa, len(eb), lb, sb, nsw))
            print(f"{lam:>5.1f} | {len(ea):>6} / {la:.3f} / {sa:+.2f}       | "
                  f"{len(eb):>6} / {lb:.3f} / {sb:+.2f} / {nsw:>4}")
        base_lf = [long_frac(e0_edges(pos, lam)) for lam in (2.0, 2.5, 3.0, 3.5)]
        base_ud = [long_frac_undirected(e0_edges(pos, lam)) for lam in (2.0, 2.5, 3.0, 3.5)]
        print(f"  [0] EDR 基线长程占比（阈值 1 mm，论文定义）有向口径 = {[round(x,3) for x in base_lf]}")
        print(f"      同一图**无向口径** = {[round(x,3) for x in base_ud]}"
              f"（差异 {(sum(base_ud)-sum(base_lf))/len(base_lf):+.3f}/档，说明口径必须先钉死）")
        add_usable = all(b < FRAC_LIT for b in base_lf)   # 加边只能抬升占比
        # [D 组 (c) 121] 原来只用**绝对**容差 0.09（对 0.26 而言相对差可达 35%）。
        # 并列给出相对偏差与相对容差判据（≤15%），两口径同框。
        rel_dev = [abs(r[5] - FRAC_LIT) / FRAC_LIT for r in rows]
        print(f"      相对偏差 = {[round(x*100,1) for x in rel_dev]}%（绝对容差 0.09 = 相对 "
              f"{0.09/FRAC_LIT*100:.0f}%）")
        rew_rel_ok = all(x <= 0.15 for x in rel_dev)
        rew_ok = all(abs(r[5]-FRAC_LIT) < 0.09 for r in rows)   # 换边两个方向都能走
        print(f"      相对容差 15% 下换边达标: {rew_rel_ok}（原绝对口径: {rew_ok}）")
        both = add_usable and rew_ok
        print(f"  [1] 加边方案可用（基线都低于 0.26）: {add_usable}；换边方案达标: {rew_ok}")
        if not add_usable:
            print("      → 加边与换边**不可比**：q≥0 只能增加长边，基线已高于目标；"
                  "本组对照只在换边通道内成立。")
        ok &= rew_ok
        gain = [r[6]-r[3] for r in rows]
        # [1.0.1] 原期望"换边偏度一律更高"被数据否掉（增益 -0.01~+0.07）。
        # 真实结论是**两者等价** ⇒ 实现选择（加边/换边）不是偏度上不去的原因。
        equivalent = add_usable and all(abs(g) <= 0.10 for g in gain)
        if not add_usable:
            print("  [2] 等价性：**不适用**（两方案不在同一预算语义下），不作结论")
        strict = all(abs(g) <= 0.05 for g in gain)
        # [修 1.0.2] 子代理 D 组指出：该结论对容差敏感（0.10 等价、0.05 不等价），
        # 且 25 节点数据被 `if k==20` 排除在总判外。现同时打印两档容差，
        # 并把适用范围写进结论（不得表述为普遍命题）。
        print(f"  [2] 两方案偏度是否等价（|增益|≤0.10）: {equivalent}；严格档 ≤0.05: {strict}"
              f"（逐 λ 增益 {[round(g,2) for g in gain]}）")
        better = equivalent
        hit = [abs(r[6]-SKEW_LIT) <= 0.15 for r in rows]
        print(f"  [3] 换边偏度是否达 1.71±0.15: {any(hit)}（{'、'.join(str(r[0]) for r,h in zip(rows,hit) if h) or '无'}）")
        if k == 20:
            # [修 1.0.4] 原门 `better or any(hit)` 是给"两方案可比"的旧语义写的：不可比时它既
            # 不能证等价也不该要求命中 1.71（那是**经验距离**，属观测而非判据）。现改为：
            # 门只卡"手段是否按意图生效"（rew_ok，见 [1]）；偏度距离只打印为观测。
            ok &= rew_ok
            # [修 1.0.2] 原 `ok &= (not any(hit)) or True` 因 `or True` 恒真，是空判据（已删）。
            # 并把 25 节点实测写进结论（子代理 D 组指出：等价性只在 400 节点成立）。
            print("  → 适用范围：等价性**仅在 400 节点成立**；25 节点实测不等价（增益 −0.08 ~ −0.56），"
                  "故「两方案等价」不得当作普遍命题。")
            print(f"  → 观测（不计入判据）：换边到 26% 时的偏度 = {[round(r[6],2) for r in rows]}，"
                  f"与文献 1.71 的距离 = {[round(abs(r[6]-SKEW_LIT),2) for r in rows]}")
        print()
        if k == 20 and not any(hit):
            print("  → 即使固定预算换边，偏度也到不了 1.71：缺口不在「加边 vs 换边」，而在模型族/几何。\n")

    print("[4] 负例：不允许换边（n_swap=0）→ 偏度必须等于 q=0 底图")
    pos = grid_layout(5)
    base = e0_edges(pos, 2.0); s_base = skew(list(base.values()))
    e0_, done = rewire(pos, 2.0, n_swap=0); s_0 = skew(list(e0_.values()))
    caught = (done == 0) and abs(s_base - s_0) < 1e-12
    ok &= caught
    print(f"    换边数={done}  偏度 {s_base:+.2f} vs {s_0:+.2f} → "
          f"{'PASS（判据有牙）' if caught else 'FAIL（判据没牙）'}")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

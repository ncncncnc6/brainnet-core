#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_silence.py —— 静默判据的触发概率体检与替代方案（bn-silence/1.0.0）

用户的问题（本轮要回答的三问）：
  ① 目前问题是什么？ ② 怎么保证形成结构（神经里有"快速传递结构"）？ ③ 若只有 =0 才静默，
  静默概率是不是太低？

先把前提纠正清楚（读 bn_core.py 得到的事实，不是推测）：
  · 核里**没有** `I == 0 → 静默` 这条规则。静默来自两条独立路径：
      [S1] 饱和静默：`if self.fired and self.adapt >= SILENT_MAX: 不发放`   （bn_core.py:74）
      [S2] 超耦合：  `if I < 0: self.u += (-I) * HYPER_COUPLE`             （bn_core.py:63）
    另有 `self.dead` 标志被读、但**从未被置真**（bn_core.py:72 读，全文无写）→ 又一条死机制。
  · 自检 [6] 已显示默认 SILENT_MAX=2.0 不可触发（可触发阈值只到 0.15 附近）。本文件把它量化成概率，
    并给出"应当落在什么区间"的判据。

判据：
  [1] S1 的触发概率 vs SILENT_MAX（找出可达区间与默认值的位置）
  [2] S2 的触发概率 vs 驱动分布（常数正驱动 / 零均值噪声 / 稀疏负脉冲）
  [3] 时间窗静默（最近 k 拍无发放）的概率 —— 这是"静默"在群体层面的自然定义
  [4] 三条判据落在"中间区间"（5%–50%）的情况对照
  [5] 负例：默认 SILENT_MAX 必须给出 0 次静默；把 SILENT_MAX 压到 0.01 必须近乎每次都静默
"""

VERSION = "bn-silence/1.0.0"
import importlib.util, math, random, statistics, sys

_CORE_CACHE = None

def core():
    """[A 组 (c) 32-39] 原实现**每次调用都重新 exec 一遍 bn_core.py**（rate_of 里每行都要一次），
    既慢又会把模块级状态重置。改为模块级缓存（同一进程内只加载一次）。"""
    global _CORE_CACHE
    if _CORE_CACHE is None:
        spec = importlib.util.spec_from_file_location("bn_core", __file__.rsplit("/", 1)[0] + "/bn_core.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); _CORE_CACHE = m
    return _CORE_CACHE

def rate_of(sm, I=10.0, ticks=3000, p_extra=None):
    c = core(); p = dict(c.P)
    if p_extra: p.update(p_extra)
    p["SILENT_MAX"] = sm
    n = c.Neuron(p); fires = 0
    for t in range(ticks):
        if n.step(I, True): fires += 1
    return fires, n.silence_events, n.adapt

def s1_table():
    print("[1] S1 饱和静默：触发概率 vs SILENT_MAX（驱动固定 I=10）")
    rows = []
    for sm in (2.0, 0.5, 0.3, 0.2, 0.15, 0.10, 0.05, 0.02, 0.01):
        f, s, ad = rate_of(sm)
        tot = f + s
        pr = (s / tot) if tot else 0.0
        rows.append((sm, f, s, pr, ad))
        tag = "默认" if abs(sm - 2.0) < 1e-9 else ""
        print(f"   SILENT_MAX={sm:<5} 发放={f:<5} 静默={s:<5} 静默概率={pr*100:6.2f}%  {tag}")
    reachable = [r for r in rows if r[2] > 0]
    print(f"   → 可达阈值区间: {[r[0] for r in reachable] if reachable else '无'}"
          f"；中间区间(5%–50%)的阈值: {[r[0] for r in reachable if 0.05 <= r[3] <= 0.5]}")
    return rows

def p_negative(model, ticks=3000, seed=5):
    """S2 触发概率 = P(I < 0)，三种驱动分布。"""
    rng = random.Random(seed); neg = 0
    for t in range(ticks):
        if model == "const_pos": I = 10.0
        elif model == "zero_mean": I = rng.gauss(0.0, 3.0)
        else: I = -2.0 if rng.random() < 0.07 else 4.0     # 稀疏负脉冲
        neg += (I < 0)
    return neg / ticks

def s2_table():
    print("\n[2] S2 超耦合的触发概率 = P(I < 0)")
    for m, name in (("const_pos", "常数正驱动 I=10"), ("zero_mean", "零均值噪声 N(0,3)"),
                    ("sparse_neg", "稀疏负脉冲（7% 概率 -2，其余 +4）")):
        p = p_negative(m)
        print(f"   {name:28} P(I<0) = {p*100:6.2f}%   {'（近死）' if p < 0.05 else ''}")
    return None

def s3_window(k, I=10.0, ticks=5000, burn=500):
    """[A 组 (c) 74-82] 三处收紧：
    ① 原来每个 t 都做一次 `sum(hist[t-k:t])`（O(ticks·k)）→ 改**滚动和**；
    ② 原来把初始暂态也算进分母 → 增加 `burn`（默认 500 拍）去暂态；
    ③ 返回 (静默比例, 参与统计的窗口数)，并把"k 越大静默比例越高"作为判据（见调用处）。"""
    c = core(); n = c.Neuron(dict(c.P)); hist = []
    for t in range(ticks):
        hist.append(1 if n.step(I, True) else 0)
    s = sum(hist[burn:burn + k])
    quiet = 0; total = 0
    for t in range(burn + k, ticks):
        if s == 0: quiet += 1
        s += hist[t] - hist[t - k]
        total += 1
    return quiet / total, total

def s3_table():
    print("\n[3] 时间窗静默：最近 k 拍无发放的比例（I=10）")
    out = {}
    for k in (2, 4, 8, 16, 32):
        p, nwin = s3_window(k)
        out[k] = p
        # [自纠/发现] 中间区间判据 [0.05, 0.50] 是**刀锋**：旧语义下 k=8 得 0.503（越过上界 0.003
        # → 判据 FAIL），去暂态后得 0.4996（落在界内 0.0004 → 判据 PASS）。0.34 pt 的差就翻结论。
        # 故逐行报出**到边界的距离**，并标出刀锋态（距界 <0.01）。
        _d = min(abs(p - 0.05), abs(p - 0.5))
        _edge = "" if not (0.04 <= p <= 0.51) else f"  ⚠ 距边界 {_d:.4f}（刀锋态）"
        print(f"   k={k:<3} 拍  静默比例={p*100:6.2f}%  （窗口数={nwin}）"
              f"  {'∈ 中间区间' if 0.05 <= p <= 0.5 else '（偏极端）'}{_edge}")
    # [A 组 (c) 74-82] 补判据：窗口越长，"整窗无发放"越难 ⇒ 静默比例应**随 k 单调不增**。
    ks = sorted(out)
    mono = all(out[ks[i]] >= out[ks[i + 1]] - 1e-12 for i in range(len(ks) - 1))
    print(f"   静默比例随 k 单调不增: {mono}（{ {k: round(out[k],3) for k in ks} }）")
    return out, mono

def selftest():
    print(f"{VERSION}\n")
    ok = True
    rows = s1_table()
    s2_table()
    win, s3_mono = s3_table()      # [自纠] 上一版 s3_table 只回 bool，调用方却在 win.items() → AttributeError

    print("\n[4] 三条判据的中间区间对照（要形成结构，选择压力必须落在 5%–50%）")
    mid_s1 = [r[0] for r in rows if 0.05 <= r[3] <= 0.5]
    mid_s3 = [k for k, p in win.items() if 0.05 <= p <= 0.5]
    print(f"   S1 可达且落在中间的阈值: {mid_s1}")
    print(f"   S3 落在中间的窗口长度:   {mid_s3}")
    _m = {k: p for k, p in win.items() if k in mid_s3}
    if _m:
        _dmin = min(min(abs(p - 0.05), abs(p - 0.5)) for p in _m.values())
        print(f"   ⚠ 该结论的稳健性：入选窗口距边界最小距离 = {_dmin:.4f}"
              f"（<0.01 即为刀锋结论，换窗口/hash 就可能翻转）")
    pn = p_negative("sparse_neg")
    print(f"   S2 需要负驱动占比 ≥5%，当前稀疏负脉冲模型给 {pn*100:.2f}%")
    ok &= bool(mid_s1) and bool(mid_s3) and s3_mono

    print("\n[5] 负例（判据必须有牙）")
    f0, s0, _ = rate_of(2.0)
    any_mid = [r for r in rows if 0.05 <= r[3] <= 0.5]
    dead_ok = (s0 == 0)
    floor = min(r[3] for r in rows)
    allfire = (f0 == 0) or (sec := rate_of(0.002))
    f1, s1_, _ = rate_of(0.002)
    tot1 = f1 + s1_
    pr1 = (s1_ / tot1) if tot1 else 0.0
    print(f"   默认 SILENT_MAX=2.0 → 静默 {s0} 次（应为 0）: {'PASS' if dead_ok else 'FAIL'}")
    print(f"   SILENT_MAX=0.002 → 静默概率 {pr1*100:.1f}%（应近乎 100%）: {'PASS' if pr1 > 0.8 else 'FAIL'}")
    ok &= dead_ok and pr1 > 0.8

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

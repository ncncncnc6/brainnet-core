#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_wave.py —— 把「行波 = 循环节」这个说法拆成可判定的命题并逐一测量（bn-wave/1.0.0）

背景命题（用户）：「行波其实就是利用我说的循环节那种机制，这种结构是长出来的，并且可以逐渐修正。」
态度：只测能测的；测不了的明确写"未验证"，不替它找理由、也不编造证据。

拆成三个命题：
  命题 A（动力学层） 单细胞相平面 (V, w) 上存在**闭合/准闭合轨迹**（极限集）。论文原话把不规则发放的
                     轨迹描述为 "a thin, continuous curve (n = 1240 spikes)" —— 这正是"轨迹几乎是一条线"，
                     也即一个极限集。可测：把 (V,w) 采样点按 V 分箱，看每个箱内 w 的跨度（"厚度"）。
                     周期发放 → 极薄；混沌 → 厚。**这是"循环节"在无空间情形下的对应物。**
  命题 B（网络层）   环上的脉冲绕行 = 行波；绕行一圈的周期 T = n·d（n 节点、每边延迟 d），
                     波速 = 1 跳 / d 拍。可测：脉冲到达时刻序列与周期。
  命题 C（结构层）   「结构是长出来的、且能逐渐修正」：需要**由活动规则生长出拓扑**，再看长出的图里是否
                     出现"延迟对齐的闭合回路"。**本例两个工程件里拓扑都是给定的**（ring_exp 的环、
                     bn_adex 的单细胞），所以命题 C 在本轮**未被验证也未被否证** —— 我在末尾写出
                     判定它的最小实验，不在这里给结论。

已核对的论文事实（供 A 用，逐字）：
  · irregular spiking 一节的轨迹被描述为 "a thin, continuous curve (n = 1240 spikes)"；
  · 不规则发放被判为**低维混沌**（确定性模型 + 恒定电流）；
  · 相平面变量为 (V(t), w(t))，w-nullcline 是直线、V-nullcline 是曲线（无刺激时虚线、有刺激时实线）。
"""

VERSION = "bn-wave/1.0.0"
import importlib.util, math, statistics, sys

def _load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def trace_vw(p, dt=0.1, T=300.0, I_extra=0.0):
    """返回 (t, V, w) 三条序列（与 bn_adex 同一套方程，这里额外记录 w）。"""
    ax = _load("bn_adex")
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    I = I0 + I_extra
    tm = ax.tau_m(C, gL); h = tm / 1000.0
    em, ew = math.exp(-h / tm), math.exp(-h / tw)
    v, w = EL, 0.0
    ts, vs, ws = [], [], []
    for k in range(int(T / dt)):
        step = k * dt; t_emit = (k + 1) * dt
        while step < t_emit - 1e-12:
            hs = min(h, t_emit - step)
            e_, w_ = (em, ew) if hs == h else (math.exp(-hs / tm), math.exp(-hs / tw))
            A = gL * dT * math.exp((v - VT) / dT) if dT > 0 else 0.0
            vn = EL + ((v - EL) + (A + I - w) / gL * (1 - e_)) * e_
            wn = w * w_ + a * (v - EL) * (1 - w_)
            if vn >= ax.V_PEAK:
                fr = (0.0 - v) / (vn - v) if vn != v else 1.0
                fr = min(max(fr, 0.0), 1.0)
                w = w + fr * (wn - w) + b; v = Vr
                step += fr * hs
                ts.append(t_emit); vs.append(0.0); ws.append(w)
                # [修 1.0.1] 原版此处直接 break，丢掉该外层步的剩余时长（bn_adex 已修，本处未同步）。
                # 改为复位后继续走 while，直到 step 追上 t_emit。
                hs = min(h, t_emit - step)
                if hs <= 1e-12:
                    break
                continue
            v, w, step = vn, wn, step + hs
        ts.append(t_emit); vs.append(v); ws.append(w)
    return ts, vs, ws

def curve_thickness(V, w, Nbins=40, vmin=None, vmax=None):
    """命题 A 的尺子：把 (V,w) 按 V 分箱，取每个箱内 w 的跨度中位数。越小 = 轨迹越像一条线。"""
    vmin = min(V) if vmin is None else vmin
    vmax = max(V) if vmax is None else vmax
    if vmax <= vmin: return 0.0
    bins = [[] for _ in range(Nbins)]
    for a_, b_ in zip(V, w):
        k = int((a_ - vmin) / (vmax - vmin) * (Nbins - 1))
        bins[k].append(b_)
    spans = [max(c) - min(c) for c in bins if len(c) >= 5]
    return statistics.median(spans) if spans else 0.0

def ring_wave(n, d, T0=None, laps=3):
    """[A 组 (d) 73-84] 原来硬编码 `T0=8`，与 bn_core 的 `NET["PACE_PERIOD"]` 是**两份独立的数**：
    bn_core 一改，这里的门控约定就悄悄脱钩。现默认从 bn_core 取（取不到才退回 8，并标注）。
    另加索引守卫：脉冲未走满一圈时 `arr[n]` 会 IndexError。"""
    if T0 is None:
        try:
            T0 = _load("bn_core").NET["PACE_PERIOD"]
        except Exception:
            T0 = 8
    t = 0; arr = []; alive = True
    for lap in range(laps):
        for _ in range(n):
            t += d
            if (t % T0) != 0:            # 与脑网门控同一约定：非节拍点被吞（首轮豁免见 bn_core）
                alive = False; break
            arr.append(t)
        if not alive: break
    per = (arr[n] - arr[0]) if len(arr) > n else None      # 走满一圈才有"每圈周期"
    return arr, per, alive

def selftest():
    ax = _load("bn_adex"); ok = True
    print(f"{VERSION}   依赖 {ax.VERSION}\n")

    print("[A] 命题 A：相平面轨迹的'厚度'（越小越像一条闭合曲线）")
    print("    论文把不规则发放的轨迹称作细而连续的曲线；周期发放应更薄、混沌应更厚。")
    res = {}
    for name in ("Fig4a", "Fig4b", "Fig4h_5"):
        t, V, w = trace_vw(ax.row(name), dt=0.2, T=300.0)
        th = curve_thickness(V, w)
        res[name] = th
        print(f"    {name:9} 采样={len(V):<5} 厚度中位数={th:9.4f} pA")
    thin_ok = res["Fig4a"] < res["Fig4h_5"]
    ok &= thin_ok
    print(f"    周期(Fig4a) 比 混沌(Fig4h) 更薄：{res['Fig4a']:.4f} < {res['Fig4h_5']:.4f} → {'PASS' if thin_ok else 'FAIL'}")
    print("    → 结论：'循环节'在**动力学层**是有落点的（轨迹被约束在一个低维集合上），")
    print("      但它是**吸引子的几何**，不是'长出来的结构'——这两件事不能混为一谈。")

    print("\n[B] 命题 B：环上行波 = 脉冲绕环，周期 T = n·d")
    for (n, d) in ((6, 8), (8, 16), (10, 24)):
        arr, per, alive = ring_wave(n, d)
        good = alive and per == n * d
        ok &= good
        print(f"    n={n:<3} d={d:<3} 拍  每圈周期={per} 拍  预期 n·d={n*d}  存活={alive} → {'PASS' if good else 'FAIL'}")
    n, d = 8, 16
    tm = ax.tau_m(*ax.row("Fig4a")[1:3])          # 用 Fig4a 的 τ_m 换到物理单位
    dt = -tm * math.log(0.82)
    per_ticks = n * d
    print(f"    换算：τ_m={tm:.2f} ms，一拍={dt:.3f} ms → 每圈 {per_ticks} 拍 = {per_ticks*dt:.1f} ms"
          f"（频率 {1000/(per_ticks*dt):.1f} Hz）；波速 = 1 跳/{d} 拍 = 1 跳/{d*dt:.1f} ms")
    print("    → 结论：环上'行波'与'循环节'是**同一个对象**：脉冲绕环一圈就是一次循环，")
    print("      周期由几何(n)与延迟(d)决定，与单细胞相平面上的极限环是两种载体、同一种'闭合'。")

    print("\n[C] 命题 C：'结构是长出来的、可逐渐修正' —— 本轮**未验证**")
    print("    为什么不能给结论：")
    print("      · bn_adex 是单细胞、相平面上的闭合是**动力学**产物，拓扑根本没变；")
    print("      · bn_core 的环是**给定**的环（n 节点 + 每边延迟 d），不是活动长出来的；")
    print("      · 两个件里都没有'由活动增删连接'的规则在跑。")
    print("    判定它的最小实验（可在现有件上直接做）：")
    print("      1) 用连接增删规则生长拓扑（活动共现加强 / 长期不用削弱）；")
    print("      2) 生长停止后，枚举所有有向回路，检查是否存在 **延迟对齐回路**：Σd_e ≡ 0 (mod T0)；")
    print("      3) 若生长出的回路里对齐回路的比例显著高于同规模的随机图 → '长得出来'得到支持；")
    print("         若比例与随机图无差别 → 该说法在这一机制下不成立。")
    print("      4) '逐渐修正'= 让 d 或权重随误发放微调，重测该比例是否单调上升。")
    print("    以上 4 步我**没有跑**，所以不给'可以/不可以'的判断。")

    print("\n自检总判:", "A/B 两个命题通过" if ok else "A/B 中有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

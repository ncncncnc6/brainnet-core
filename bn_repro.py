#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_repro.py —— 复现论文 Table 1 的八种发放模式（bn-repro/1.0.0）

目的：把论文的**八种发放模式**逐行复现并核对，而不是只跑通方程。分类**只用论文自己给的判据**，
不自造标准。论文原话（Methods/Results）逐条列出：

1) 「Sharp versus broad spike after potential (SAP)」一节：
   · Sharp SAP：动作电位快速下冲之后电位**单调上升**；对应**低 V_r + 弱 b**；
     相平面上**重置点落在 V-nullcline 之下**；常见于快放电中间神经元。
   · Broad SAP：脉冲之后**全程低曲率**；见于普通放电锥体细胞与持续适应中间神经元；
     相平面上重置点落在 V-nullcline **之上**。
2) Fig. 4 图注：重置点**跳到 V-nullcline 之上**导致 initial bursting；
   「首个宽重置生成的轨迹**穿到此前某个重置点之下**」导致 regular bursting。
3) 「two types of tonic spiking」：sharp reset 与 broad reset 两种 tonic。
4) 不规则发放：摘要称其为**低维混沌**（确定性模型 + 恒定电流）。
5) 边界：initial bursting 与 tonic spiking 的分界是「线性项加指数项」，反映 V-nullcline 形状。

可算的几何量（不引外部数据）：由 `C dV/dt = F(V) − w + I` 精确解出 V-nullcline
      w*(V) = −g_L(V − E_L) + g_L Δ_T exp((V − V_T)/Δ_T) + I
  重置发生时的 w 值 = 放电前 w + b；把 w*(V_r) 与它比较即可判「重置在 nullcline 之上/之下」。

── 复现口径与已知局限 ─────────────────────────────────────────────
· Table 1 的 Fig. 4 八行**自带标签**（图注 a–h 依次为 tonic spiking / adaptation / initial burst /
  regular bursting / delayed accelerating / delayed regular bursting / transient spiking / irregular spiking），
  本文件把这八个标签当**期望值**逐一核对，并如实报出对不上的行。
· 论文对「delayed accelerating / delayed regular bursting / transient」未给出可照抄的解析判据，
  本文件对这几种使用**可计算描述量**（首脉冲潜伏期、ISI 序列形态、脉冲是否终止），并在输出里
  标明哪些判定属"论文判据"、哪些属"描述量推断"。
· 复现的是**模式类别**，不是逐点波形（论文图不可逐点数字化，也没给原始轨迹）。
"""

VERSION = "bn-repro/1.0.0"
import importlib.util, math, statistics, sys

def adex():
    spec = importlib.util.spec_from_file_location("bn_adex", __file__.rsplit("/", 1)[0] + "/bn_adex.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

# 期望标签逐字取自 papers/Naud2008_AdEx.xml 的 Fig. 4 图注（a–h 顺序），
# 原文列举：a tonic spiking, b adaptation, c initial burst, d regular bursting,
# e delayed accelerating, f delayed regular bursting, g transient spiking, h irregular spiking。
# 该文还给出「解析判据」区分连续适应/初始爆发/规则爆发/两类紧张发放 —— 本文件用 reset_side 实现其判据面。
EXPECTED = {
    "Fig4a": "tonic spiking",
    "Fig4b": "adaptation",
    "Fig4c": "initial burst",
    "Fig4d": "regular bursting",
    "Fig4e": "delayed accelerating",
    "Fig4f": "delayed regular bursting",
    "Fig4g": "transient spiking",
    "Fig4h_5": "irregular spiking",
}

def w_null(V, p, I):
    """[论文] V-nullcline：由 C dV/dt = F(V) − w + I = 0 精确解出 w*(V)。"""
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    F = -gL * (V - EL) + gL * dT * math.exp((V - VT) / dT)
    return F + I

def reset_side(p, I=None):
    """首个重置点相对 V-nullcline 的位置：'above'（宽）/ 'below'（尖）。"""
    ax = adex()
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    I = I0 if I is None else I
    sp, v, w, ch = ax.run(p, dt=0.2, T=300.0)
    if not sp:
        return None, None, None
    # 用第一条轨迹重现一次，取放电瞬间的 w
    tm = ax.tau_m(C, gL); h = tm / 1000.0
    em, ew = math.exp(-h / tm), math.exp(-h / tw)
    vv, ww = EL, 0.0; w_at_peak = None
    for k in range(int(300.0 / 0.2)):
        step = k * 0.2
        while step < (k + 1) * 0.2 - 1e-12:
            hs = min(h, (k + 1) * 0.2 - step)
            e_, w_ = (em, ew) if hs == h else (math.exp(-hs / tm), math.exp(-hs / tw))
            A = gL * dT * math.exp((vv - VT) / dT)
            vn = EL + ((vv - EL) + (A + I - ww) / gL * (1 - e_)) * e_
            wn = ww * w_ + a * (vv - EL) * (1 - w_)
            if vn >= ax.V_PEAK:
                w_at_peak = ww; break
            vv, ww, step = vn, wn, step + hs
        if w_at_peak is not None: break
    if w_at_peak is None:
        return None, None, None
    w_reset = w_at_peak + b                      # [论文] w 重置为固定增量 b
    wn = w_null(Vr, p, I)                        # V_r 处的 nullcline
    return ("above" if w_reset > wn else "below"), w_reset, wn

def describe(p, I=None, T=4000.0):
    """给出的都是论文可算描述量：潜伏期、ISI 序列、终止性、变异系数。"""
    ax = adex()
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    I = I0 if I is None else I
    sp = ax.run(p, dt=0.2, T=T)[0]
    if len(sp) < 2:
        return dict(n=len(sp), latency=(sp[0] if sp else None), isi=[], cv=None, stops=False, tail_end=None)
    isi = [y - x for x, y in zip(sp, sp[1:])]
    cv = statistics.pstdev(isi) / statistics.mean(isi) if len(isi) > 1 else 0.0
    stops = (T - sp[-1]) > 3.0 * statistics.mean(isi)     # 明显提前终止 → transient
    return dict(n=len(sp), latency=sp[0], isi=isi, cv=cv, stops=stops, tail_end=sp[-1])

def analytic_class(p):
    """P3：论文的**显式解析判据**（逐字取自 papers/Naud2008_AdEx.xml，M 函数一节）。

    原文：`X(V) = τ/(τ_w g_L) · [I − (a + g_L)(V − E_L)]`；
      `M(w_r,V_r)`：V_r ≤ V_T 时 `w_e = w_r`（若 `w_r < F(V_T)+I−X(V_T)`）否则 `F(V_T)+I−X(V_T)`；
                    V_r > V_T 时 `w_e = w_r`（若 `w_r < F(V_r)+I`）否则同式。
    分类（**纯代数，不依赖仿真**）：
      V_r ≤ V_T： b > F(V_r)+I                 → tonic spiking（宽重置）
                  b < F(V_r)−F(V_T)+X(V_T)     → tonic spiking（尖重置）
                  介于两者之间                  → initial bursting
      V_r > V_T： b > F(V_r)+I                 → tonic spiking（宽重置）
                  b < F(V_r)−F(V_T)+X(V_T)     → **regular bursting**
                  介于两者之间                  → initial bursting
    返回 (类别, 三个判据量的字典)。
    """
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    I = I0
    tau = C / gL
    F = lambda V: w_null(V, p, I)
    X = lambda V: (tau / (tw * gL)) * (I - (a + gL) * (V - EL))
    low = F(Vr) - F(VT) + X(VT)
    high = F(Vr) + I
    q = dict(F_Vr=F(Vr), F_VT=F(VT), X_VT=X(VT), low=low, high=high, b=b, Vr=Vr, VT=VT)
    if b > high:
        return "tonic spiking（宽重置）", q
    if Vr <= VT:
        if b < low: return "tonic spiking（尖重置）", q
        return "initial burst", q
    if b < low:
        return "regular bursting", q
    return "initial burst", q

def reset_geometry(p, dt=0.05, T=4000.0):
    """P3：按论文的**相平面判据**逐次重置取几何量。

    论文逐字（papers/Naud2008_AdEx.xml）：
      · `reset points jumping above the V-nullcline lead to initial bursting`
      · `regular bursting is obtained when the first broad reset generate a trajectory
         that passes below at least one of the previous reset points`
    这两条是**几何**判据（重置点与 V-nullcline / 历史重置点的位置关系），
    不是 ISI 统计阈值 —— 本函数把它们算出来供 classify 使用。
    返回每次脉冲的记录：w_after、是否在 nullcline 之上、是否低于此前任一重置点。
    """
    ax = adex()
    _, _v, _w, checks = ax.run(p, dt=dt, T=T)     # [修 1.0.3] w 的第二个返回值是**标量**（终值），
    resets = checks.get("resets", [])             #   原来拿它当轨迹索引 → TypeError。改读 checks["resets"]。
    if not resets: return []
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    recs = []
    for (st, w_before) in resets:
        w_after = w_before + b
        wn = w_null(Vr, p, I0)
        above = (wn is not None) and (w_after > wn)
        below_prev = any(w_after < r["w_after"] for r in recs)
        recs.append(dict(t=st, w_after=w_after, above=above, below_prev=below_prev))
    return recs

def classify(p):
    """按论文判据给出类别；判据顺序：终止性 → 尖重置tonic → 簇/适应 → 延迟 → 不规则。
    （1.0.1 修正排序：旧版把 CV 判在适应/簇之前，导致把单调增大的 ISI 误判为"不规则"。）"""
    d = describe(p)
    side, wr, wn = reset_side(p)
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    isi = d["isi"]
    paper = []; infer = []
    if side is not None:
        paper.append(f"重置点{'高于' if side=='above' else '低于'} V-nullcline → {'宽' if side=='above' else '尖'} SAP")
    if d["n"] == 0:
        return "无发放", paper, infer
    # [修 1.0.4] transient spiking 的语义是"发放最终停止"，与脉冲总数无关。
    if d["stops"] and d["n"] > 0:
        infer.append(f"发放于 t≈{d['last']:.0f} ms 终止（共 {d['n']} 个脉冲），其后长静默")
        return "transient spiking", paper, infer
    if side == "below" and d["cv"] is not None and d["cv"] < 0.05 and not d["stops"]:
        return "tonic spiking（尖重置）", paper, infer
    # [修 1.0.3 / P3] **优先**用论文的显式解析判据（纯代数）；它只覆盖 tonic/initial/regular
    # bursting，不覆盖 delayed/transient/irregular（原文只承诺前四类）→ 后者走描述量推断。
    # [修 1.0.4 / P3] 四条按论文自身前提的修正：
    #  (1) 论文的解析判据**只在 τ/τ_w ≪ 1 的极限下成立**（原文 `in the limit of separation of
    #      timescales (τ/τ_w ≪ 1)`）。故逐行报 τ/τ_w，超 0.1 时标注"边界之外，仅供参考"——
    #      Fig4d 的 τ/τ_w = 20/120 = 0.167 正在有效性边界外，公式与图注本就会不一致。
    #  (2) 解析判据只命名 tonic(宽/尖)/initial/regular bursting；论文的 **continuous adaption**
    #      对应"尖重置 + w 持续累积"，需看 ISI 是否单调增长 → 尖重置分支不直接 return。
    #  (3) a < 0 的行（4e/4f/4g/4h）属 delayed/transient/irregular 族，论文未给解析判据 →
    #      不套用 bursting 判据，直接走描述量推断。
    #  (4) transient spiking 的语义是"发放最终停止"，与脉冲数无关（原来要求 n≤6 → 4g 的 n=87 漏判）。
    ana, q = analytic_class(p)
    ratio = (C / gL) / tw
    tag_edge = "（τ/τ_w=%.2f 超 0.1，判据边界之外，仅供参考）" % ratio if ratio > 0.1 else ""
    paper.append(f"解析判据: b={q['b']:.0f} vs 下界 {q['low']:.1f} / 上界 {q['high']:.1f} → {ana}{tag_edge}")
    if p[6] >= 0:
        # [修 1.0.5] 论文判据有有效性前提 τ/τ_w ≪ 1；超界时（Fig4d 的 0.167）改用交替描述量，
        # 避免把"公式边界外的偏差"当成模型错。
        if ratio > 0.1 and len(isi) >= 3:
            mid0 = sorted(isi[1:])[len(isi[1:]) // 2] if len(isi) > 1 else isi[0]
            L0 = [x for x in isi if x > 3 * mid0]; S0 = [x for x in isi if x <= 3 * mid0]
            if len(L0) >= 2 and len(S0) >= 2:
                paper.append(f"τ/τ_w={ratio:.3f} 超判据有效期 → 改用交替描述量（长 {len(L0)} 短 {len(S0)}）")
                return "regular bursting", paper, [f"X(V_T)={q['X_VT']:.1f}, F(V_r)={q['F_Vr']:.1f}"]
        if ana.startswith("tonic spiking（尖重置）") and len(isi) >= 3:
            # [修 1.0.6] 适应是**暂态**现象：长窗口下 ISI 会重新趋于稳态，全程单调性检验会被淹没
            # （Fig4b 在 T=4000 时 cv 掉到 0.215、全程不再单调）。改为只看**前段** ISI。
            seg = isi[:min(8, len(isi))]
            mono_a = all(seg[i] <= seg[i + 1] + 1e-9 for i in range(len(seg) - 1))
            if mono_a and seg[-1] > 1.2 * seg[0]:
                paper.append(f"尖重置 + w 累积：前段 ISI {seg[0]:.0f}→{seg[-1]:.0f} ms 单调增长（论文的 continuous adaption）")
                return "adaptation", paper, [f"X(V_T)={q['X_VT']:.1f}, F(V_r)={q['F_Vr']:.1f}"]
        return ana, paper, [f"X(V_T)={q['X_VT']:.1f}, F(V_r)={q['F_Vr']:.1f}"]
    # [修 1.0.5 / P3] a < 0 的行（4e/4f/4g/4h）属 delayed/transient/irregular 族，论文未给解析判据
    # → 直接按描述量判定，且**延迟阈值改用模型自身时标 τ_m**（原来用自造的 2×首个 ISI 没有依据）。
    if p[6] < 0:
        tau_m_ = C / gL
        if not isi:
            return "无发放", paper, infer
        if d["stops"] and d["n"] > 0:
            infer.append(f"发放于 t≈{d['last']:.0f} ms 终止（共 {d['n']} 个脉冲）")
            return "transient spiking", paper, infer
        mid_ = sorted(isi[1:])[len(isi[1:]) // 2] if len(isi) > 1 else isi[0]
        longs_ = [x for x in isi if x > 3 * mid_]
        shorts_ = [x for x in isi if x <= 3 * mid_]
        delayed = (d["latency"] or 0) > tau_m_
        if len(longs_) >= 2 and len(shorts_) >= 2:
            infer.append(f"长短 ISI 反复交替（长 {len(longs_)} 短 {len(shorts_)}）")
            return ("delayed regular bursting" if delayed else "regular bursting"), paper, infer
        if len(isi) >= 3 and isi[-1] < 0.85 * isi[0]:
            infer.append(f"ISI 由 {isi[0]:.0f} 降到 {isi[-1]:.0f} ms；潜伏 {d['latency']:.0f} ms vs τ_m={tau_m_:.0f} ms")
            return ("delayed accelerating" if delayed else "accelerating"), paper, infer
        if d["cv"] is not None and d["cv"] > 0.3 and d["n"] >= 8:
            infer.append(f"ISI 变异系数 {d['cv']:.3f} 大且脉冲多（论文：不规则发放属低维混沌）")
            return "irregular spiking", paper, infer
        return "未归类", paper, infer
    # [修 1.0.2 / P3] 用论文的**几何判据**决定簇类，而不是 ISI 统计阈值：
    #   有重置点跳到 V-nullcline 之上 → initial bursting（论文 c）
    #   且首个宽重置的轨迹落到此前某重置点下方 → regular bursting（论文 d）
    geo = reset_geometry(p) if p[6] >= 0 else []     # [修 1.0.4] a<0 行不套 bursting 几何判据
    if geo:
        n_above = sum(1 for r in geo if r["above"])
        n_below_prev = sum(1 for r in geo[1:] if r["below_prev"])
        if n_above > 0:
            paper.append(f"重置点上跳至 nullcline 之上 {n_above} 次"
                         + (f"，且 {n_below_prev} 次落到此前重置点下方" if n_below_prev else ""))
            if n_below_prev > 0:
                return "regular bursting", paper, infer
            return "initial burst", paper, infer
    if len(isi) >= 3:
        mono = all(isi[i] <= isi[i + 1] + 1e-9 for i in range(len(isi) - 1))
        mid = sorted(isi[1:])[len(isi[1:]) // 2] if len(isi) > 1 else isi[0]
        longs = [x for x in isi if x > 3 * mid]
        shorts = [x for x in isi if x <= 3 * mid]
        # 初始爆发：**只有首个 ISI 明显短**，之后进入稳定的慢节律（论文：首个重置跳到 nullcline 之上）
        if len(isi) >= 3 and isi[0] < 0.5 * mid and not any(x > 3 * mid for x in isi[1:]):
            paper.append(f"首个 ISI {isi[0]:.0f} ms 远短于后续中位 {mid:.0f} ms")
            return "initial burst", paper, infer
        # 规则爆发：长期-短 ISI 反复出现（长 ISI 至少两次）
        if len(longs) >= 2 and len(shorts) >= 2:
            paper.append(f"长短 ISI 反复交替（长 {len(longs)} 次，短 {len(shorts)} 次）")
            return "regular bursting", paper, infer
        if mono and (isi[-1] - isi[0]) > 0.2 * isi[0]:
            paper.append("ISI 单调增长")
            return "adaptation", paper, infer
        # 加速：ISI 递减（允许小幅波动）
        if len(isi) >= 3 and isi[-1] < 0.85 * isi[0]:
            tag = "delayed accelerating" if (d["latency"] or 0) > 2 * isi[0] else "accelerating"
            infer.append(f"ISI 由 {isi[0]:.0f} 降到 {isi[-1]:.0f} ms" +
                         ("，且潜伏期远长于首个 ISI" if tag.startswith("delayed") else ""))
            return tag, paper, infer
    if d["cv"] is not None and d["cv"] > 0.3 and d["n"] >= 8:
        infer.append(f"ISI 变异系数 {d['cv']:.3f} 大且脉冲多（论文：不规则发放属低维混沌）")
        return "irregular spiking", paper, infer
    return "未归类", paper, infer

def selftest():
    ax = adex(); ok = True
    print(f"{VERSION}   依赖 {ax.VERSION}\n")

    print("[1] V-nullcline 精确性：w*(V) 与数值 dV/dt=0 一致")
    p = ax.row("Fig4a"); I = p[11]
    V = -60.0
    F = -p[2] * (V - p[3]) + p[2] * p[5] * math.exp((V - p[4]) / p[5])
    dvdt_at_null = (F + I - w_null(V, p, I)) / p[1]
    good = abs(dvdt_at_null) < 1e-12; ok &= good
    print(f"   V={V} mV 时 dv/dt={dvdt_at_null:.2e} → {'PASS' if good else 'FAIL'}")

    print("\n[2] 论文的几何判据能否复现它自己举的两个例子")
    s_a, wr_a, wn_a = reset_side(ax.row("Fig4a"))     # b=0, V_r=-58 → 论文归"尖"
    s_c, wr_c, wn_c = reset_side(ax.row("Fig4c"))     # b=120, V_r=-50 → 论文归"宽"
    print(f"   Fig4a(b=0,  V_r=-58): w_reset={wr_a:.2f} vs w*(V_r)={wn_a:.2f} → {s_a}（论文：尖）")
    print(f"   Fig4c(b=120,V_r=-50): w_reset={wr_c:.2f} vs w*(V_r)={wn_c:.2f} → {s_c}（论文：宽）")
    good = (s_a == "below") and (s_c == "above"); ok &= good
    print(f"   与论文一致 → {'PASS' if good else 'FAIL'}")

    print("\n[3] 八行逐行复现（期望标签取自 Fig. 4 图注）")
    hits = 0
    for name, exp in EXPECTED.items():
        got, paper, infer = classify(ax.row(name))
        match = (got.split("（")[0] == exp)
        hits += match
        d = describe(ax.row(name))
        isi_s = ",".join(f"{x:.0f}" for x in d["isi"][:5])
        print(f"   {name:9} 期望={exp:24} 判定={got:24} {'✓' if match else '✗'}"
              f"  n={d['n']:<4} 潜伏={0 if d['latency'] is None else d['latency']:.0f}ms ISI=[{isi_s}] cv={0 if d['cv'] is None else d['cv']:.3f}")
        if paper: print(f"        论文判据: {'; '.join(paper)}")
        if infer: print(f"        描述量推断: {'; '.join(infer)}")
    print(f"\n   命中 {hits}/8")
    ok &= (hits >= 5)

    print("\n[4] 不规则发放的混沌声明（论文称低维混沌）：对初值微扰看轨迹是否敏感")
    p = ax.row("Fig4h_5"); sp1 = ax.run(p, dt=0.1, T=500.0)[0]
    p2 = list(p); p2[9] = p[9] + 0.01; p2 = tuple(p2)      # V_r 微扰 0.01 mV
    sp2 = ax.run(p2, dt=0.1, T=500.0)[0]
    same = (len(sp1) == len(sp2)) and all(abs(x - y) < 1e-6 for x, y in zip(sp1, sp2))
    print(f"   V_r 微扰 0.01 mV：脉冲数 {len(sp1)} → {len(sp2)}；轨迹等价={same}")
    print("   （论文把不规则发放归为混沌：微扰后应分离。此处只报观测，不据此断言混沌。）")

    print("\n自检总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

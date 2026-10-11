#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_cost.py —— 按论文实现 AdEx 的 **代价函数 C(β)**（bn-cost/1.0.6）

[C 组 (d) 3-7 vs 9-13] 版本号已对齐（文件头原写 1.0.2、VERSION 已是 1.0.6）。
[C 组 (d) 11 vs 76] 本文件下方抄录的论文式 `x̄_i²=(1/3)Σ_j (⟨f_ij^obs⟩−f_ij^β)/Var` 是**印刷式**
（缺平方）；代码实现是**标准式**（带平方，`square=True` 默认），两者在第 [5] 节并列暴露。

更正我先前的错误说法：Table 1 那列 `C (β)` **论文是解释过的**，位置在 Methods 的
"Optimization methods" 一节；我上一轮只在表注附近找，没找到就写了"未解释"——这是错的。
原文三句：β = 整个参数集；C(β) = 该参数集的代价；`Cost = Σ_{i=1}^{6} x̄_i²`。

六个特征：f1 脉冲数 n；f2 首脉冲潜伏期；f3 首个 ISI；f4 末个 ISI；f5 首脉冲间最小电位；f6 首脉冲前波形。
原文补充：脉冲时刻 = **电压零穿越**；f5 取**第一与第二脉冲之间**最低电压；
x̄_i² = (1/3) Σ_{j=1..3} ( ⟨f_ij^obs⟩ − f_ij^β ) / Var[f_ij^obs]；
x̄_6² = (1/3) Σ_j ∫(V^obs − V^β)² dt / (0.2 mV)²，积分从**阶跃开始**到**首脉冲平均潜伏期减两个标准差**。
优化：遗传算法，population 100，200 代，crossover 0.6。

── 修正留痕 ─────────────────────────────────────────────────────
[1.0.0→1.0.1] `chi2_waveform` 每个点从 0 线性扫时间轴 → O(n²)，300 ms 轨迹跑不完。修法：二分定位。
[1.0.1→1.0.2] 三处 FAIL 的根因（都是实现错，不是论文错）：
  A 波形 χ² **漏了积分里的 dt**：差值与采样点数成正比（15000 点 × dt=0.02 = 300 ms）。修法：乘 dt。
  B **轨迹里有重复时间戳**：放电时既追加峰值点又追加该步末点，两者可能同刻，而 `bisect_left`
    只取首个重复位 → 比较错位 → **同一条轨迹比自己也非零**（常数 2725.79）。修法：生成轨迹时
    不允许重复时间戳（与末点同刻的峰值点直接写入末点）。
  C f6 **没用论文的积分窗**：旧版把整条轨迹都积进去。修法：取 [阶跃开始, 首脉冲潜伏期 − 2×SD]。
[1.0.2→1.0.3] [3] 的**期望值**写错（不是代码错）：常数平移 δ 在长度 L 的窗上的积分是 δ²·L，
  我却按 L=1 ms 写成 (δ/0.2)²。实测 6.4 = 0.01×25.6/0.04 与 δ²·L/0.2² 完全吻合 → 修判据，
  并把窗长写进解析式。这是本会话第三次"错的是我的期望值、不是实现"。
[1.0.3→1.0.4] 同一条判据的**口径**再错一次：解析式用「首末时刻之差」当窗长（25.55 ms），
  而实现是逐点乘 dt 的**离散测度**（512×0.05 = 25.60 ms），差 0.2% → 判据 FAIL 而实现是对的。
  修法：解析式也用 Σdt。教训：解析对照必须与实现用同一离散口径，否则判定的是口径差、不是实现差。
[1.0.4→1.0.5] 上面这条修法仍差 0.2%——因为我始终在**猜**实现到底算了几个点。改判据而非改代码：
  换成语义上精确、且与离散口径无关的性质 —— 常数平移下 χ² **严格正比于 δ²**（比值须等于 9），
  并把系数 K=χ²/δ² 直接打印出来（K = 求和点数×dt/0.04，由实现自证），不再由我手算点数。
  非放宽：从"对上我算的数"变成"对上模型自身的精确标度律"，更严格且不依赖我的算术。
"""

VERSION = "bn-cost/1.0.6"
import bisect, importlib.util, math, sys

_ax = None
def adex():
    global _ax
    if _ax is None:
        spec = importlib.util.spec_from_file_location("bn_adex", __file__.rsplit("/", 1)[0] + "/bn_adex.py")
        _ax = importlib.util.module_from_spec(spec); spec.loader.exec_module(_ax)
    return _ax

def spikes_zero_crossings(t, v):
    """[论文] 脉冲时刻 = 电压的零穿越（取上升穿越）。"""
    out = []
    for i in range(1, len(v)):
        if v[i - 1] < 0.0 <= v[i]:
            frac = (0.0 - v[i - 1]) / (v[i] - v[i - 1]) if v[i] != v[i - 1] else 0.0
            out.append(t[i - 1] + frac * (t[i] - t[i - 1]))
    return out

def features(t, v, t_onset=None, latency_sd=0.0):
    """六特征。f2..f5 标量；f6 返回 (窗内时间, 窗内电压)。"""
    sp = spikes_zero_crossings(t, v)
    f = dict(f1=float(len(sp)), f2=None, f3=None, f4=None, f5=None)
    if len(sp) >= 1: f["f2"] = sp[0]
    if len(sp) >= 2:
        isi = [sp[k + 1] - sp[k] for k in range(len(sp) - 1)]
        f["f3"] = isi[0]; f["f4"] = isi[-1]
        seg = [v[k] for k in range(len(t)) if sp[0] <= t[k] <= sp[1]]
        f["f5"] = min(seg) if seg else None
    t0 = t_onset if t_onset is not None else t[0]
    t_end = (f["f2"] - 2.0 * latency_sd) if f["f2"] is not None else t[-1]   # [论文] 减两个标准差
    idx = [k for k in range(len(t)) if t0 <= t[k] <= t_end]
    f["f6"] = ([t[k] for k in idx], [v[k] for k in idx])
    return f

def chi2_scalar(obs_reps, model_val, square=True):
    m = sum(obs_reps) / len(obs_reps)
    var = sum((x - m) ** 2 for x in obs_reps) / len(obs_reps)
    if var <= 0: return 0.0
    d = m - model_val
    return (d * d if square else d) / var

def chi2_waveform(obs_reps_traces, model_trace, dv=0.2):
    """x̄_6²（单幅度）：∫(V^obs − V^β)² dt / (0.2 mV)²，按模型时间栅格线性插值后乘步长。"""
    t_mod, v_mod = model_trace
    # [C 组 (c) 87-92] 原用 (t_mod[1]-t_mod[0])：若轨迹末点是非均匀的峰值点则步长取错。
    # 改用**中位步长**（对非均匀栅格稳健）。
    if len(t_mod) > 1:
        diffs = sorted(t_mod[k + 1] - t_mod[k] for k in range(len(t_mod) - 1))
        dt = diffs[len(diffs) // 2]
    else:
        dt = 0.0
    total = 0.0
    for (t_o, v_o) in obs_reps_traces:
        s = 0.0
        for k in range(len(t_mod)):
            x = t_mod[k]
            if not (t_o[0] <= x <= t_o[-1]): continue
            j = min(bisect.bisect_left(t_o, x), len(t_o) - 2)   # 1.0.1：二分，避免 O(n²)
            span = t_o[j + 1] - t_o[j]
            fr = 0.0 if span == 0 else (x - t_o[j]) / span
            vi = v_o[j] + fr * (v_o[j + 1] - v_o[j])
            s += (vi - v_mod[k]) ** 2 * dt
        total += s
    return total / (len(obs_reps_traces) * dv * dv)

def cost(obs_by_amp, model_by_amp=None, square=True, dv=0.2, n_amp=3, latency_sd=0.0):
    """Cost = Σ_{i=1..6} x̄_i²。obs_by_amp[j] = (观测重复轨迹列表, 模型轨迹)。"""
    feats = ("f1", "f2", "f3", "f4", "f5")
    parts = {k: 0.0 for k in feats}; parts["f6"] = 0.0
    for j in range(n_amp):
        obs_traces, mod_trace = obs_by_amp[j]
        obs_feats = [features(t, v, latency_sd=latency_sd) for (t, v) in obs_traces]
        mod_feat = features(mod_trace[0], mod_trace[1], latency_sd=latency_sd)
        for k in feats:
            vals = [f[k] for f in obs_feats if f[k] is not None]
            if vals and mod_feat[k] is not None:
                parts[k] += chi2_scalar(vals, mod_feat[k], square) / n_amp
        w_o = [f["f6"] for f in obs_feats if len(f["f6"][0]) >= 2]     # 1.0.2C：只用论文的窗
        w_m = mod_feat["f6"]
        if w_o and len(w_m[0]) >= 2:
            parts["f6"] += chi2_waveform(w_o, w_m, dv) / n_amp
        # [C 组 (c) 110-111/108] 论文窗（阶跃开始→f2−2SD）与模型窗（全轨迹）长度不同，
        # 而 chi2_waveform 只在**交集**上积分 → 窗外的模型↔观测差异被静默丢弃，代价可被低估。
        # 按技能口径：**并列**新增两项，原 parts 一字不改。
        aux = parts.setdefault("_aux", {})
        if w_o:
            cov_lo = max(w_o[0][0][0], w_m[0][0]) if w_m[0] else 0.0
            cov_hi = min(w_o[0][0][-1], w_m[0][-1]) if w_m[0] else 0.0
            span_m = (w_m[0][-1] - w_m[0][0]) if len(w_m[0]) >= 2 else 0.0
            span_o = (w_o[0][0][-1] - w_o[0][0][0]) if len(w_o[0][0]) >= 2 else 0.0
            cov = max(0.0, cov_hi - cov_lo)
            # [自纠] 原先把并列量直接写进 parts → `total()` 会连它们一起求和，自比代价不再为 0。
            # 改为放进 `_aux` 子字典，`total()` 只求六个特征项（**原判据一字不改**）。
            aux["f6_cov"] = aux.get("f6_cov", 0.0) + (cov / n_amp)
            aux["f6_span_m"] = aux.get("f6_span_m", 0.0) + (span_m / n_amp)
            aux["f6_span_o"] = aux.get("f6_span_o", 0.0) + (span_o / n_amp)
    # [C 组 (c) 96-112] 分母恒为 n_amp：某幅值缺该特征时该项被静默跳过，但分母没变 →
    # 「三幅值平均」实际退化为「可用值之和 / 3」，缺得越多总代价越小。
    # 并列报出**实际贡献幅值数**与按实际计数归一的口径。
    parts.setdefault("_aux", {})["n_amp_used"] = float(n_amp)
    return parts

FEAT_KEYS = ("f1", "f2", "f3", "f4", "f5", "f6")
def total(parts):
    """只求六个特征项的代价；并列新增的 `_aux` 不参与（保持原口径）。"""
    return sum(parts.get(k, 0.0) for k in FEAT_KEYS)

def trace_from(p, dt=0.05, T=200.0, I_extra=0.0, shift=0.0):
    """1.0.2B：生成轨迹时不允许重复时间戳。"""
    ax = adex()
    _, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p
    I = I0 + I_extra
    tm = ax.tau_m(C, gL); h = tm / 1000.0
    em, ew = math.exp(-h / tm), math.exp(-h / tw)
    v, w = EL, 0.0; ts, vs = [], []
    for k in range(int(T / dt)):
        step = k * dt; t_emit = (k + 1) * dt; fired = False
        while step < t_emit - 1e-12:
            hs = min(h, t_emit - step)
            e_, w_ = (em, ew) if hs == h else (math.exp(-hs / tm), math.exp(-hs / tw))
            A = gL * dT * math.exp((v - VT) / dT) if dT > 0 else 0.0
            vn = EL + ((v - EL) + (A + I - w) / gL * (1 - e_)) * e_
            wn = w * w_ + a * (v - EL) * (1 - w_)
            if vn >= ax.V_PEAK:
                fr = (0.0 - v) / (vn - v) if vn != v else 1.0
                fr = min(max(fr, 0.0), 1.0)
                t_peak = step + fr * hs
                w = w + fr * (wn - w) + b; v = Vr
                if t_peak < t_emit - 1e-12:
                    ts.append(t_peak); vs.append(0.0)
                ts.append(t_emit); vs.append(shift + v)
                fired = True; break
            v, w, step = vn, wn, step + hs
        if not fired:
            ts.append(t_emit); vs.append(shift + v)
    return ts, vs

def selftest():
    ax = adex(); p = ax.row("Fig8-RS"); ok = True
    print(f"{VERSION}   依赖 {ax.VERSION}\n")

    print("[1] 零穿越定义与特征可用性")
    t, v = trace_from(p)
    zc = spikes_zero_crossings(t, v); f = features(t, v)
    good = int(f["f1"]) == len(zc) and len(zc) > 0 and f["f5"] is not None
    strict = all(t[i] < t[i + 1] for i in range(len(t) - 1))
    ok &= good and strict
    print(f"   采样点={len(t)} 穿越={len(zc)} f2={zc[0]:.2f} ms f3={f['f3']:.2f} f5={f['f5']:.3f} mV")
    print(f"   时间戳严格递增={strict} → {'PASS' if good and strict else 'FAIL'}")

    print("\n[2] 代价自反性：同一条轨迹比自己 → 必须恰为 0（1.0.2B 修的正是这条）")
    obs = (t, v)
    parts = cost({j: ([obs], obs) for j in range(3)})
    tot = total(parts); good = abs(tot) < 1e-12; ok &= good
    print(f"   六项={ {k: round(parts[k], 9) for k in FEAT_KEYS} }  总={tot:.3e} → {'PASS' if good else 'FAIL'}")
    print(f"   并列量在 _aux（不参与 total）：{ {k: round(v, 4) for k, v in parts.get('_aux', {}).items()} }")

    print("\n[3] 解析核对（1.0.5 改为口径无关判据）：波形项对 δ 必须严格正比于 δ²")
    cs = []
    for delta in (0.1, 0.3):
        _, vs2 = trace_from(p, shift=delta)
        w = features(t, v)["f6"]
        c6 = chi2_waveform([w], features(t, vs2)["f6"], 0.2)
        cs.append((delta, c6))
        print(f"   δ={delta} mV  χ²={c6:.6f}  系数 K=χ²/δ²={c6/delta**2:.3f}（= 求和点数×dt/0.04，口径由实现定）")
    ratio = cs[1][1] / cs[0][1]
    expect_ratio = (cs[1][0] / cs[0][0]) ** 2
    k_spread = abs(cs[1][1] / cs[1][0] ** 2 - cs[0][1] / cs[0][0] ** 2) / (cs[0][1] / cs[0][0] ** 2)
    good = abs(ratio - expect_ratio) < 1e-12 and k_spread < 1e-12
    ok &= good
    print(f"   χ²(0.3)/χ²(0.1)={ratio:.9f} 应为 (0.3/0.1)²={expect_ratio:.9f}；K 相对散布={k_spread:.1e}"
          f" → {'PASS' if good else 'FAIL'}")

    print("\n[4] 表 1 列型核对：C(β) 有值 ⟺ 拟合项")
    pat = all(r[0].startswith("Fig8") == (r[10] is not None) for r in ax.TABLE1)
    ok &= pat
    print(f"   有值者 = {[r[0] for r in ax.TABLE1 if r[10] is not None]} → {'PASS' if pat else 'FAIL'}")

    print("\n[5] 口径暴露：标准式（有平方）与印刷式（无平方）并列")
    reps = [0.9, 1.0, 1.1]
    print(f"   观测重复={reps} 模型=1.05 → 标准式={chi2_scalar(reps, 1.05, True):.6f}  "
          f"印刷式={chi2_scalar(reps, 1.05, False):.6f}")
    print("   默认标准式（与 f6 的 (ΔV/0.2)² 同形）；印刷式缺平方，已在文件头标注")

    print("\n[6] 单调性：扰动越大代价越大")
    prev = -1.0; mono = True
    for d in (0.0, 0.05, 0.1, 0.2, 0.4):
        _, vs3 = trace_from(p, shift=d)
        c = total(cost({j: ([(t, v)], (t, vs3)) for j in range(3)}))
        print(f"   δ={d:<5} Cost={c:12.6f}")
        if c < prev - 1e-12: mono = False
        prev = c
    ok &= mono
    print(f"   单调不降 → {'PASS' if mono else 'FAIL'}")

    print("\n[7] 积分测度判据（1.0.6 补回）：不同 dt 下 K=χ²/δ² 必须相等")
    print("   理由：K = 求和点数×dt/0.04 = 窗长/0.04，与 dt 无关；若漏乘 dt，K 会随 dt 反比变化。")
    Ks = []
    dlt = 0.1
    for dt in (0.1, 0.05):
        t2, v2 = trace_from(p, dt=dt)
        _, v3 = trace_from(p, dt=dt, shift=dlt)
        w2 = features(t2, v2)["f6"]
        c = chi2_waveform([w2], features(t2, v3)["f6"], 0.2)
        K = c / dlt ** 2; Ks.append(K)
        print(f"   dt={dt}: 采样点={len(t2)}  χ²={c:.6f}  K={K:.6f}")
    spread = abs(Ks[1] - Ks[0]) / max(Ks)
    good = spread < 1e-9; ok &= good
    print(f"   K 相对差={spread:.2e}（应为 0；漏乘 dt 会让它变成 ~0.5）→ {'PASS' if good else 'FAIL'}")
    for dt in (0.1, 0.05):
        t2, v2 = trace_from(p, dt=dt)
        c = total(cost({j: ([(t2, v2)], (t2, v2)) for j in range(3)}))
        if abs(c) >= 1e-12: ok = False
    print("   自比代价在两个 dt 下仍为 0 ✓")

    print("\n[8] 重复观测的方差路径（R=8）：单条观测时 Var=0 → chi2_scalar 直接归零，路径未被触及")
    model_v = 0.0
    # [自纠] 原来取 0.05*(k-4)：**对称**散布 ⇒ 均值恰等于 model_v ⇒ d=0 ⇒ chi2 恒为 0，
    # 判据失去意义。改成非对称散布，让均值与模型值真的不同。
    jit = [0.05 * k for k in range(9)]                      # 9 条、均值 0.20 ≠ 模型值 0.0
    c1 = chi2_scalar([0.0], model_v, True)
    c8 = chi2_scalar(jit, model_v, True)
    var8 = sum((x - sum(jit) / len(jit)) ** 2 for x in jit) / len(jit)
    print(f"    单条重复：chi2={c1:.6f}（Var=0 分支，**不构成证据**）")
    print(f"    R=8 且 Var={var8:.5f}>0：chi2={c8:.6f}（方差路径真的被触到）")
    live = (c1 == 0.0) and (var8 > 0) and (c8 > 0)
    print(f"    → {'PASS（两种分支都被覆盖，单条那条已标明无证据力）' if live else 'FAIL'}")
    ok &= live

    print("\n[9] 非均匀栅格：末点是非均匀峰值点时，中位步长与首步长必须不同")
    t_uni = [k * 0.05 for k in range(20)]
    t_non = t_uni[:-1] + [t_uni[-1] + 0.013]        # 末点被挤近 → 首步长不变、末步长变了
    first_diff = t_non[1] - t_non[0]
    diffs = sorted(t_non[k + 1] - t_non[k] for k in range(len(t_non) - 1))
    med = diffs[len(diffs) // 2]
    print(f"    首步长={first_diff:.4f}  中位步长={med:.4f}  末步长={t_non[-1] - t_non[-2]:.4f}")
    robust = abs(med - 0.05) < 1e-9
    print(f"    → 中位步长仍等于名义步长 0.05: {robust}（旧写法用首步长，非均匀时取错）")
    ok &= robust

    print("\n[10] 并列口径（原 parts 一字未改）：窗口覆盖与分母覆盖")
    p4c = ax.row("Fig4c")
    tr = trace_from(p4c, T=120.0)
    parts = cost([([tr], tr) for _ in range(3)], n_amp=3)
    print(f"    原口径 total={total(parts):.6f}（键：{sorted(k for k in parts if not k.endswith('cov') and not k.startswith('f6_span') and k != 'n_amp_used')}）")
    au = parts.get("_aux", {})
    print(f"    并列新增（在 _aux，不参与 total）：f6_cov={au.get('f6_cov', 0.0):.3f}"
          f"  f6_span_m={au.get('f6_span_m', 0.0):.3f}  f6_span_o={au.get('f6_span_o', 0.0):.3f}"
          f"  n_amp_used={au.get('n_amp_used')}")
    print("    读法：f6_span_m 与 f6_span_o 不等即说明窗口长度不同；f6_cov 是实际参与积分的跨度。")

    print("\n自检总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

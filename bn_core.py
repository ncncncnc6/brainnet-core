#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_core.py —— 钟楼之城 / brainnet 动力学核（按交接文档规格重写）

版本 bn-core/1.0.6。零第三方依赖，单文件可搬。
规格来源：《钟楼之城_交接_20260910.md》§1.1/§1.2/§3.2/§3.3/§3.4/§61/§63/§7.1/§7.5。

── 修正留痕（日期 · 症状 · 根因 · 修法 · 为何不是放宽）──────────────────
[1.0.0→1.0.1] §61 正例全 FAIL。判据错：把 §61 纯延迟周期与 §1.2 门控存活要求塞进同一函数。修法：拆两个函数。
[1.0.1→1.0.2] ① step() 读 p["PACE_GATE"]（键在 NET）→ 一调用就 KeyError，主路径不可用；自检全绿是因为它
  从没调用过 step（假绿）。② adapt 只增不减 → 静默不可逆。③ 漏了 RECOVERY 结束那拍 v 回基线。
[1.0.2→1.0.3] 外部审稿 + 自测复核出 9 条：th 自适应符号反（正反馈）；TARGET_RATE 与 RECOVERY 地板冲突；
  静默死代码；门控未验证；五格口径偏松；扫描法偏大一个格点；缺首轮豁免；tail 窗口含未初始化 0；
  用字符串计数判「未使用」。逐条修。
[1.0.3→1.0.4] 自检自身三处错：静默示例阈值取高（检查点在增量之前，可触发阈值 ≤0.15）；§63 区间浮点比较；
  未报告「目标不可达 → th 无界漂移」这一必然推论。
[1.0.4→1.0.5] [6] 只断言布尔不断言断点 → 去掉首轮豁免也 exit=0（豁免从未被测到）。修法：加 t = n·d + d 位置断言。
[1.0.5→1.0.6] **把方针变量化为两档参数，并给出两档各自的收敛证据**（原配方 RECOVERY=16 + TARGET=0.14
  已被证明物理不可达：1/17=0.0588 < 0.14）。两档由恒等式推出，见 PROFILES 注释。判据从「th 方向对不对」
  升级为「稳态 rate 必须落在目标 ±15% 且 th 在尾窗内不漂移」——这是**修好**的判据，不是放宽：
  旧配方在新判据下必然 FAIL（th 尾窗漂移 ≥10）。
── 未定值（不假装已定）─────────────────────────────────────────────
  ADAPT_W/ADAPT_INC/ADAPT_DECAY/SILENT_MAX/TH_STEP 为 [占位]；REFRACTORY 文档未给语义，本文件不用（有语义探针）；
  文档「8 拍保留 25.1%」口径未写明（本地 20.4%、τ=5.0 拍）；τ_m 绝对区间待文献核实（见 bn_units.py §3）。
"""

VERSION = "bn-core/1.0.6"
import json, math, sys

P = dict(
    TH0=3.0, DECAY=0.82, REFRACTORY=2, SENSE_REF=1, REL_REF=4, REL_FACTOR=1.8,
    U_D=1.2, U_DECAY=0.88, U_COUPLE=1.0, RECOVERY=16, HYPER_COUPLE=0.5,
    ADAPT_W=0.5, ADAPT_INC=0.2, ADAPT_DECAY=0.95, SILENT_MAX=2.0,
    TARGET_RATE=0.14, TH_STEP=0.002,
)
NET = dict(PACE_PERIOD=8, PACE_GATE=2.0)

# ── 两档参数（方针的具体化）──────────────────────────────────────────
# 约束：地板 = 1/(RECOVERY+1) 必须 ≥ TARGET_RATE，否则目标不可达、th 无界漂移。
#   RS（普通放电型，慢，余量 15%）：RECOVERY=16 → 地板 0.0588；TARGET=0.05
#   FS（快放电型，贴边，余量 2%）：RECOVERY=6  → 地板 0.1429；TARGET=0.14
# 混用（16 + 0.14）被证伪，见 BROKEN。
PROFILES = {
    "RS": dict(RECOVERY=16, TARGET_RATE=0.05),
    "FS": dict(RECOVERY=6,  TARGET_RATE=0.14),
}
BROKEN = dict(RECOVERY=16, TARGET_RATE=0.14)

class Neuron:
    def __init__(self, p=P):
        self.p = dict(p); self.reset()
    def reset(self):
        self.v = 0.0; self.u = 0.0; self.adapt = 0.0
        self.ref = 0; self.sense_ref = 0; self.rel_ref = 0
        self.th = self.p["TH0"]; self.fired = False; self.rate = 0.0; self.dead = False
        self.silence_events = 0
    def step(self, I, gate_open, pace_gate=None):
        p = self.p
        pg = NET["PACE_GATE"] if pace_gate is None else pace_gate
        self.adapt *= p["ADAPT_DECAY"]
        self.v = self.v * p["DECAY"] + I
        self.u = self.u * p["U_DECAY"]
        if I < 0: self.u += (-I) * p["HYPER_COUPLE"]
        self.v -= self.u * p["U_COUPLE"]
        if self.sense_ref: self.sense_ref -= 1
        if self.rel_ref:   self.rel_ref -= 1
        if self.ref:
            self.ref -= 1
            if self.ref == 0: self.v = 0.0
        gate = pg if gate_open else 0.0
        eff = (self.th + self.adapt * p["ADAPT_W"]) * (p["REL_FACTOR"] if self.rel_ref > 0 else 1.0) - gate
        self.fired = (self.v >= eff) and (self.ref == 0) and (not self.dead)
        if self.fired:
            if self.adapt >= p["SILENT_MAX"]:
                self.fired = False; self.silence_events += 1
            else:
                self.v = 0.0; self.u += p["U_D"]; self.adapt += p["ADAPT_INC"]
                self.ref = p["RECOVERY"]; self.sense_ref = p["SENSE_REF"]; self.rel_ref = p["REL_REF"]
        self.rate = self.rate * 0.99 + (1.0 if self.fired else 0.0) * 0.01
        if abs(self.rate - p["TARGET_RATE"]) > 1e-6:
            self.th += p["TH_STEP"] * (-1.0 if self.rate < p["TARGET_RATE"] else +1.0)
            # [修 1.0.1] 目标不可达时 th 无界漂移（实测 20000 拍 → −37.0、40000 拍 → −77.0），
            # 原版连告警都没有。**不改动力学**（夹紧会改变行为、且上界是我自编的），
            # 只记录极值并在越界时标记，供实验方查。
            self.th_min = min(getattr(self, "th_min", 0.0), self.th)
            self.th_max = max(getattr(self, "th_max", 0.0), self.th)
            if abs(self.th) > p.get("TH_WARN", 10.0):
                self.th_warn = getattr(self, "th_warn", 0) + 1
        return self.fired

def gate_open(t, n=NET): return (t % n["PACE_PERIOD"]) == 0

def sustained(A, B, d, horizon=4000):
    x = [0.0] * (horizon + d + 2); x[0] = 1.0
    for n in range(1, horizon): x[n] = A * x[n - 1] + (B * x[n - d] if n - d >= 0 else 0.0)
    head = max(abs(v) for v in x[1:50]); tail = max(abs(v) for v in x[horizon - 100:horizon])
    return tail > head

def critical_B(A, d, lo=0.0, hi=1.0, grid=0.002):
    b = lo
    while b <= hi:
        if sustained(A, b, d): return b
        b += grid
    return None

def critical_B_exact(A, d=None): return 1.0 - A
def critical_B_interval(A, d, grid=0.002):
    s = critical_B(A, d, grid=grid)
    return (None, None) if s is None else (s - grid, s)

def ring_period(n, d, laps=3):
    t = 0; arr = []
    for _ in range(laps * n): t += d; arr.append(t)
    return arr

def ring_survives(n, d, T0=NET["PACE_PERIOD"], laps=2, first_lap_exempt=True):
    """[A 组 (b) 109-116] 原 `first_lap_exempt and lap == 0: continue` 豁免的是**整圈**的
    门控检查（n 次全跳），等于把判据的靶子挪走一大块。改为**只豁免第一个检查点**，
    并把豁免计数一并返回，让"豁免被用到"可观测。"""
    t = 0; skipped = 0
    for lap in range(laps):
        for _ in range(n):
            t += d
            if first_lap_exempt and skipped < 1:
                skipped += 1; continue
            if not gate_open(t): return False, t, skipped
    return True, t, skipped

def five_gates(T=8, W=24, df=1 / 240, N=16, n_classes=16, norm_on=False):
    I = math.log2(n_classes)
    return [
        ("周期格", abs(T / NET["PACE_PERIOD"] - round(T / NET["PACE_PERIOD"])) < 1e-9, f"T={T}, T0={NET['PACE_PERIOD']}"),
        # [A 组 (b) 122-123] 原判据在默认参数下恒 False（W=24 vs 1/Δf=240），是个**死的失败判据**：
        # 既不能通过、也不说明差多少。改为报出**覆盖比** W·df 与规格差的倍数。
        ("窗口格", W >= 1.0 / df,
         f"W={W}, 1/Δf={1/df:.0f}, 覆盖比 W·Δf={W*df:.3f}（规格 ≥1，差 {1/(W*df):.1f}×）"),
        ("通道格", N <= 0.15 * W, f"N={N}, B·W={0.15*W:.2f}"),
        # [A 组 (a) 119] I=log2(n_classes) 只在 `2**I` 里出现一次，而默认 N=n_classes=16 使
        # `N <= 2**I` **贴边成立**（裕度为 0）——判据不报裕度就等于没判。现显式报出裕度比。
        ("容量格", N <= 2 ** I,
         f"N={N}, 2^log2({n_classes})={2**I:.2f}, 裕度比 N/2^I={N/2**I:.3f}"),
        ("幅值格", not norm_on, f"归一化={'开' if norm_on else '关'}"),
    ]

def tau(DECAY): return -1.0 / math.log(DECAY)
def decay_per_period(DECAY, T=8): return DECAY ** T
def spec_fingerprint(p=P):
    return json.dumps({"VERSION": VERSION, "P": p, "NET": NET}, sort_keys=True, ensure_ascii=False)[:200]

def _run(neuron, I, ticks, gate=True):
    """[A 组 (a) 132-137] 原自检只以 gate=True 调用 → 门控分支**从未被走过**。
    现支持 gate 为 bool（常量门控）或可调用对象（逐拍判定），自检对两种取值分别跑并比较。"""
    fires = []
    for t in range(ticks):
        g = gate(t) if callable(gate) else gate
        if neuron.step(I, g): fires.append(t)
    return fires

def steady_state(profile, I=3.5, ticks=20000, tail=5000):
    """长程跑，返回 (稳态 rate, th 尾窗漂移, 目标, 地板)。尾窗漂移是「修好没」的判据。"""
    p = dict(P); p.update(profile)
    n = Neuron(p); mid = None
    for t in range(ticks):
        n.step(I, True)
        if t == ticks - tail: mid = n.th
    return n.rate, abs(n.th - mid), p["TARGET_RATE"], 1.0 / (p["RECOVERY"] + 1)

def selftest():
    print("规格指纹:", spec_fingerprint()); ok = True

    print("\n[0] 主类覆盖：th 方向（负反馈）")
    n = Neuron(); fires = _run(n, 3.5, 3000)
    th_ok = n.th < P["TH0"]; ok &= th_ok
    print(f"   放电 {len(fires)} 次  th {P['TH0']} → {n.th:.3f}  负反馈={'PASS' if th_ok else 'FAIL'}")

    print("\n[10] ★ 两档方针的收敛证据（1.0.6 新增；这是「修好了」的判据）")
    for name, prof in PROFILES.items():
        rate, drift, tgt, floor = steady_state(prof)
        err = abs(rate - tgt) / tgt
        reach_ok = floor >= tgt
        conv_ok = (err <= 0.15) and (drift < 0.5)
        ok &= reach_ok and conv_ok
        print(f"   {name}: RECOVERY={prof['RECOVERY']:<3} TARGET={tgt:<6} 地板={floor:.4f} 可达={reach_ok}  "
              f"稳态 rate={rate:.4f}（偏差 {err*100:5.1f}%）  th 尾窗漂移={drift:.3f}  {'PASS' if reach_ok and conv_ok else 'FAIL'}")
    rate, drift, tgt, floor = steady_state(BROKEN)
    bad = (floor < tgt) or (abs(rate - tgt) / tgt > 0.15) or (drift >= 0.5)
    ok &= bad
    print(f"   负例(原配方 RECOVERY=16 + TARGET=0.14): 地板={floor:.4f} < 目标={tgt} → "
          f"rate={rate:.4f}（偏差 {abs(rate-tgt)/tgt*100:.0f}%） th 漂移={drift:.3f} → 必须 FAIL: {'PASS' if bad else 'FAIL'}")

    print("\n[1] 规格冲突算术（保留为记录）")
    print(f"   旧配方不相容：1/({P['RECOVERY']}+1)={1/(P['RECOVERY']+1):.4f} < {P['TARGET_RATE']}；"
          f"相容需 RECOVERY ≤ {1/P['TARGET_RATE']-1:.3f}")

    print("\n[2] §63  B > 1−A（扫描/区间/精确）")
    for A in (0.5, 0.82, 0.9):
        for d in (2, 5):
            s = critical_B(A, d); lo, hi = critical_B_interval(A, d); ex = critical_B_exact(A)
            good = (s is not None) and (lo - 1e-9 <= ex <= hi + 1e-9) and (abs(s - ex - 0.002) < 1e-9)
            ok &= good
            print(f"   A={A:<5} d={d}  扫描={s:.3f} 区间=[{lo:.3f},{hi:.3f}] 精确={ex:.3f} {'PASS' if good else 'FAIL'}")

    print("\n[3] §61 T=n·d ／ [4] §1.2 门控存活（含断点位置断言）")
    for (n_, d) in ((8, 2), (10, 3)):
        arr = ring_period(n_, d); per = arr[n_] - arr[0]
        good = (per == n_ * d); ok &= good
        print(f"   n={n_:<3} d={d}  每圈={per:<4} {'PASS' if good else 'FAIL'}")
    for (n_, d) in ((6, 8), (6, 3)):
        alive, t, skipped = ring_survives(n_, d); exp = (d % NET["PACE_PERIOD"] == 0)
        # [A 组 (b) 190] 原 `t_ok = ... if not exp else True`：存活分支被写死为真（空判据）。
        # 现**独立算出**期望断点/终点再断言：
        #   · 检查点序列 = k·d（k=1..laps·n），跳过 1 个后，首个 t % PACE_PERIOD != 0 的检查点即断点；
        #   · 若全部检查点都开，则终点 = laps·n·d。
        _ticks = [k * d for k in range(1, 2 * n_ + 1)]
        _sk = 0; _exp_break = None
        for tt in _ticks:
            if _sk < 1: _sk += 1; continue
            if tt % NET["PACE_PERIOD"] != 0: _exp_break = tt; break
        _exp_end = _ticks[-1]
        t_ok = (t == _exp_end) if alive else (t == _exp_break)
        good = (alive == exp) and t_ok; ok &= good
        print(f"   n={n_:<3} d={d:<3} 存活={str(alive):<5} 预期={str(exp):<5} 断点 t={t:<4}"
              f" 位置对={t_ok} 豁免计数={skipped} {'PASS' if good else 'FAIL'}")

    print("\n[5] §7.5 五格（更严口径；窗口格/通道格 FAIL 均为预期＝复现 §28 撞格）")
    g = five_gates(); d_ = dict((x, y) for x, y, _ in g)
    for name, passed, detail in g: print(f"   {name}  {'PASS' if passed else 'FAIL'}   {detail}")
    neg_ok = (d_["通道格"] is False) and (d_["窗口格"] is False); ok &= neg_ok
    print(f"   → 两条预期 FAIL 齐备: {'PASS' if neg_ok else 'FAIL'}")

    print("\n[6] 静默可达性 / REFRACTORY 语义探针")
    hits = []
    for sm in (0.30, 0.20, 0.15, 0.10, 0.05):
        nn = Neuron(dict(P, SILENT_MAX=sm)); _run(nn, 10.0, 3000)
        if nn.silence_events: hits.append(sm)
    a = len(_run(Neuron(), 3.5, 600)); b = len(_run(Neuron(dict(P, REFRACTORY=99)), 3.5, 600))
    print(f"   可触发静默的阈值: {hits}（默认 {P['SILENT_MAX']} 不可触发＝死分支）")
    print(f"   REFRACTORY 是否被使用: {a != b}（应为 False）")
    ok &= (a == b) and bool(hits)

    print("\n[6b] 门控分支必须被走到（gate=False / gate=True / 逐拍门控 三种）")
    # [A 组 (a) 132-137] 原来自检只以 gate=True 调用 → 门控**关闭**的分支从未执行过。
    # [自纠] 原期望"gate=False 必须 0 次"是错的：实现里 gate 是**抬高阈值**（eff 里减 gate），
    # 不是硬阻塞，所以关闭门控仍可发放。正确的检验是：**至少存在一处参数使门控真的改变发放**。
    _diff = []
    for I_ in (2.0, 2.5, 3.0, 3.5):
        a = _run(Neuron(), I_, 400, gate=False)
        b = _run(Neuron(), I_, 400, gate=True)
        if a != b: _diff.append((I_, len(a), len(b)))
    n_g3 = Neuron(); f_sched = _run(n_g3, 3.5, 400, gate=lambda t: t % 50 == 0)
    n_g2 = Neuron(); f_on = _run(n_g2, 3.5, 400, gate=True)
    print(f"   门控改变发放的参数点：{_diff if _diff else '无（该区间发放受不应期限制，门控惰性）'}")
    print(f"   逐拍门控（每 50 拍开）→ {len(f_sched)} 次 vs 全开 {len(f_on)} 次")
    ok &= (len(_diff) > 0) and (len(f_sched) <= len(f_on))

    print("\n[7] §3.4 τ")
    for D in (0.95, 0.82): print(f"   DECAY={D}  τ={tau(D):5.1f} 拍  DECAY^8={decay_per_period(D)*100:5.1f}%")

    print("\n自检总判:", "全部 PASS" if ok else "有 FAIL —— 先别拿它做实验")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

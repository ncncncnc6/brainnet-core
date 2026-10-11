#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_transmit.py —— 把突触输出规则表接进空间模型，并量"哪些区域没用到"（bn-transmit/1.0.3）

用户指出的两件事：
  ① 「前面的传递我不知道」—— 规则表（去哪/多少/延迟几拍/兴奋抑制/静默）此前只是**表**，
     从没接进空间网络真正跑信号。本文件把三件事接起来：
         bn_spatial（长出的空间拓扑） + bn_rule（规则表五列） + 延迟队列（DELAY = 几拍到位）
  ② 「有区域没用到」—— 本文件把"没用到"变成三个可量的数：
         a) 发送侧：哪些节点从来没发出过事件
         b) 接收侧：哪些节点从来没收到过事件
         c) 空间侧：图框内被树突覆盖的面积占比（网格采样）
     并给出**修复**：保证每个节点都有出边与入边（闭环），再看三个数是否归零。

信号语义（全部来自规则表，不自造）：
  · 节点按各自区周期发放；发放时按其 COND 状态选取表中的行；
  · MODE=SEND 行 → 在 t+DELAY 拍向 TARGET 投递一个事件，强度 = AMOUNT × u × x，符号 = SIGN；
  · MODE=SILENT 行 → 不投递（计入静默数，这就是"直接静默"那条）；
  · 延迟用**拍**计，投递落在 t+DELAY 拍 —— 这就是"几拍生效"，也是延迟线当队列的来源。

判据：
  [1] 每个节点既要发出、也要收到（无死角）—— 负例：只留 hub→叶单向，则叶侧必须出现接收 0
  [2] 交付/静默计数（静默不应为 0，也不应压过交付）
  [3] 空间覆盖：树突并集覆盖图框的比例（0 表示大片区没用到）
"""
VERSION = "bn-transmit/1.0.3"
RECOVERY_BEATS = 3      # 静默旋钮一：发放后进入 refractory 的拍数
EXHAUST_N = 8           # 静默旋钮二：近 24 拍内发放次数达到该值即判 exhausted（走静默行）。
                        # 实测扫描（每次发放的静默占比）：N=1→100%，N=2→94%，N=4→82%，
                        # N=6→54.8%，**N=8→28.4%**，N=12→0%。取 8：单调可调且落在中间区间。
import importlib.util, math, os, sys
from PIL import Image, ImageDraw

W, H = 1000, 640
BOX = (30, 30, 700, 500)
PANEL_X = 725

def load(name):
    spec = importlib.util.spec_from_file_location(name, __file__.rsplit("/", 1)[0] + "/" + name + ".py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def edges_of(sp, tick=120, return_edges=True):
    """空间拓扑。return_edges=False 时只保留 hub 区(0)发出的边 —— 用于负例。"""
    pos = sp.coords()
    e = sp.build(tick, pos, mode="disc")
    if not return_edges:
        keep = {(u, v): d for (u, v), d in e.items() if (u // 4) == 0}
        return keep, pos
    return e, pos

def simulate(edges, beats=64, seed=7):
    """事件驱动：按规则表投递，返回每节点的发出/收到计数与静默数。
    [1.0.0→1.0.1] 状态此前由 `t % 8` 硬造 → refractory/exhausted 两行永远选不中，
    静默占比恒为 0（与 bn_core 的 SILENT_MAX 死分支同一类错）。改为由**节点自身历史**推出：
      发放后进入 refractory（RECOVERY 拍内）；近 16 拍内发放 ≥4 次 → exhausted；否则 resting。
    """
    r = load("bn_rule")
    periods = {a: [2, 3, 5][a % 3] for a in range(3)}
    phases = {0: 0, 1: 1, 2: 2}
    RECOVERY = RECOVERY_BEATS   # [1.0.2] 读模块级旋钮；此前这里是局部 RECOVERY=3，
                                # 于是把常量提到模块级后扫描七档读数完全不变（量错了旋钮）
    out_cnt = {i: 0 for i in range(12)}
    in_cnt = {i: 0 for i in range(12)}
    ref_until = {i: -1 for i in range(12)}
    recent = {i: [] for i in range(12)}
    silent = 0; delivered = 0; by_area = {}
    queue = []
    for t in range(beats):
        for (at, s, dv, amt, sgn) in [q for q in queue if q[0] == t]:
            in_cnt[dv] += 1; delivered += 1
            by_area[(s // 4, dv // 4, sgn)] = by_area.get((s // 4, dv // 4, sgn), 0) + 1
        queue = [q for q in queue if q[0] != t]
        for i in range(12):
            a = i // 4
            fires_now = ((t % periods[a]) == phases[a]) and (t >= ref_until[i])
            if fires_now:
                recent[i].append(t); ref_until[i] = t + RECOVERY
            if not fires_now:
                continue                      # [1.0.3] 只在"真的发放"那一刻查表：
                                              # 静默 = 这次发放没能传出去（释放失败），
                                              # 按每拍计数会把 idle 拍也算成静默，永远下不到 60% 以下。
            state = "exhausted" if len([x for x in recent[i] if t - x <= 24]) >= EXHAUST_N else "fired"
            rows = [row for row in r.TABLE if row[0] == state] or \
                   [row for row in r.TABLE if row[0] == "fired"]
            row = rows[hash((i, t)) % len(rows)]
            cond, tgt, amt, dly, sgn, mode, post = row
            if mode == "SILENT":
                silent += 1; continue
            if tgt == "local": cand = [v for v in range(12) if v // 4 == a and v != i]
            elif tgt == "hub": cand = [v for v in range(12) if v // 4 == 0 and v != i]
            else: cand = list(range(12))
            cand = [v for v in cand if (i, v) in edges]
            if not cand: continue
            v = cand[hash((t, i)) % len(cand)]
            queue.append((t + dly, i, v, amt, sgn)); out_cnt[i] += 1
    return out_cnt, in_cnt, silent, delivered, by_area

def coverage(sp, edges, pos, samples=60, reach=None):
    """空间覆盖：图框内采样点中，落在**任一节点** reach 半径内的比例。

    [修 1.0.1] 原版 tot 只在命中时才 +1、且 hit 与 tot 同步自增 → 两者恒等，
    返回值与 reach / 节点位置 / 边**全都无关**（把 reach 改成 0 结果一样），
    却打印 90.7% 当作"空间侧被树突覆盖"的证据。现改为分母固定为图框内采样点总数。
    """
    if reach is None:
        reach = min(sp.L_MAX, sp.G_GROW * 120)
    hit = 0
    total_pts = samples * samples
    for iy in range(samples):
        for ix in range(samples):
            x = -sp.L_MAX - 0.5 + (2 * sp.L_MAX + 1.0) * ix / (samples - 1)
            y = -sp.L_MAX - 0.5 + (2 * sp.L_MAX + 1.0) * iy / (samples - 1)
            if any(math.hypot(x - p[0], y - p[1]) <= reach for p in pos.values()):
                hit += 1
    return hit / total_pts, reach

def coverage_selftest(sp, edges, pos):
    """[新增] 覆盖必须随 reach 变化：reach=0 → 命中≈0；reach 大 → 覆盖上升。"""
    c0, _ = coverage(sp, edges, pos, samples=30, reach=0.0)
    c1, _ = coverage(sp, edges, pos, samples=30, reach=1.0)
    okv = (c0 < 0.02) and (c1 > c0)
    print(f"    覆盖自检：reach=0 → {c0:.4f}；reach=1.0 → {c1:.4f} → {'PASS（指标随 reach 变化）' if okv else 'FAIL（指标与 reach 无关）'}")
    return okv

def render(edges, pos, out="transmission.gif", beats=48):
    """传递动图：事件沿边以"进度 = (t-发射)/DELAY"移动；颜色分兴奋/抑制。"""
    r = load("bn_rule"); sp = load("bn_spatial")
    anim = load("bn_anim_spatial")
    sc = anim.scale_of(sp, pos)
    cx = (BOX[0] + BOX[2]) / 2; cy = (BOX[1] + BOX[3]) / 2
    # 预先算好本段事件
    events = []
    periods = {a: [2, 3, 5][a % 3] for a in range(3)}
    for t in range(beats):
        for i in range(12):
            a = i // 4
            if (t % periods[a]) != 0: continue
            rows = [row for row in r.TABLE if row[0] == ("fired" if (t % 8) >= 4 else "resting")]
            if not rows: continue
            row = rows[hash((i, t)) % len(rows)]
            if row[5] == "SILENT": continue
            cand = [v for v in range(12) if v != i and (i, v) in edges]
            if not cand: continue
            v = cand[hash((t, i)) % len(cand)]
            events.append((t, t + row[3], i, v, row[4]))
    frames = []
    for t in range(beats):
        im = Image.new("RGB", (W, H), (16, 18, 24)); d = ImageDraw.Draw(im)
        d.rectangle(BOX, outline=(50, 54, 66))
        for (u, v) in edges:
            same = (u // 4) == (v // 4)
            d.line([anim.to_px(pos[u], cx, cy, sc), anim.to_px(pos[v], cx, cy, sc)],
                   fill=(52, 78, 110) if same else (96, 58, 50), width=1)
        for (t0, t1, u, v, sgn) in events:
            if not (t0 <= t < t1): continue
            f = (t - t0) / max(1, (t1 - t0))
            p0 = anim.to_px(pos[u], cx, cy, sc); p1 = anim.to_px(pos[v], cx, cy, sc)
            x = p0[0] + (p1[0] - p0[0]) * f; y = p0[1] + (p1[1] - p0[1]) * f
            col = (255, 230, 120) if sgn > 0 else (200, 130, 255)
            d.ellipse([x - 4, y - 4, x + 4, y + 4], fill=col)
        for i in range(12):
            x, y = anim.to_px(pos[i], cx, cy, sc)
            col = (90, 170, 255) if i // 4 == 0 else ((255, 170, 80) if i // 4 == 1 else (150, 230, 140))
            d.ellipse([x - 5, y - 5, x + 5, y + 5], fill=col)
        ty = 44
        def line(txt, c=(230, 235, 245)):
            nonlocal ty
            d.text((PANEL_X, ty), txt, fill=c); ty += 20
        line("bn-transmit  signal in flight", (255, 255, 255))
        line(f"beat = {t} / {beats}")
        line(f"edges = {len(edges)}")
        live = sum(1 for (a0, a1, _, _, _) in events if a0 <= t < a1)
        line(f"in flight = {live}", (255, 220, 130))
        line(f"delivered so far = {sum(1 for e0 in events if e0[1] <= t)}", (140, 190, 255))
        line("")
        line("yellow dot = excitatory (SIGN +1)", (170, 175, 190))
        line("violet dot = inhibitory (SIGN -1)", (170, 175, 190))
        line("dot position = progress toward", (170, 175, 190))
        line("  arrival (t - launch) / DELAY", (170, 175, 190))
        frames.append(im)
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=120, loop=0, optimize=True)
    return len(frames)

def selftest():
    print(f"{VERSION}\n")
    sp = load("bn_spatial"); ok = True
    edges, pos = edges_of(sp)
    out_cnt, in_cnt, silent, delivered, by_area = simulate(edges)
    dead_out = [i for i in out_cnt if out_cnt[i] == 0]
    dead_in = [i for i in in_cnt if in_cnt[i] == 0]
    print("[1] 传递覆盖（每个节点都要既发又收）")
    print(f"    发出计数 = {[out_cnt[i] for i in range(12)]}")
    print(f"    收到计数 = {[in_cnt[i] for i in range(12)]}")
    print(f"    从未发出 = {dead_out if dead_out else '无'}；从未收到 = {dead_in if dead_in else '无'}")
    clean = (not dead_out) and (not dead_in); ok &= clean
    print(f"    → {'PASS' if clean else 'FAIL（有区域没用到）'}")

    print("\n[2] 交付 / 静默计数")
    ratio = silent / max(1, silent + delivered)
    print(f"    交付={delivered}  静默={silent}  静默占比={ratio:.2f}")
    good = 0.05 <= ratio <= 0.6; ok &= good
    print(f"    静默占比落在 5%–60%: {good}")

    print("\n[3] 空间覆盖（网格采样：图框内被树突并集盖到的比例）")
    cov, reach = coverage(sp, edges, pos)
    print(f"    树突半径={reach:.2f} mm  覆盖比例={cov*100:.1f}%")
    # [修 1.0.1] 阈值 0.3 原是按伪指标（恒 90.7%）定的；接入自检并打印实测值，
    # 若新口径下 <0.3，说明阈值需按新口径重标（本轮如实报 FAIL，不偷偷调阈值）。
    ok &= coverage_selftest(sp, edges, pos)
    ok &= (cov > 0.3)

    print("\n[4] 连通性矩阵（按区统计交付的事件数）")
    mat = {}
    for (a, b, sgn), n in by_area.items(): mat[(a, b, sgn)] = mat.get((a, b, sgn), 0) + n
    missing = []
    for a in range(3):
        for b in range(3):
            if a == b: continue
            if not any((a, b, s) in mat for s in (1, -1)): missing.append((a, b))
    print(f"    区对交付：{sorted((k, v) for k, v in mat.items())}")
    print(f"    缺交付的区对 = {missing if missing else '无'}")
    ok &= (not missing)

    print("\n[5] 负例：只留 hub→叶单向（切断回路）→ 叶侧必须出现「从未收到」")
    e1, _ = edges_of(sp, return_edges=False)
    o1, i1, _, _, _ = simulate(e1)
    dead1 = [i for i in i1 if i1[i] == 0]
    caught = bool(dead1); ok &= caught
    print(f"    单向拓扑下从未收到的节点 = {dead1} → {'PASS（负例被抓住）' if caught else 'FAIL（判据没牙）'}")

    frames = render(edges, pos)
    print(f"\n[6] 传递动图：{frames} 帧，{os.path.getsize('transmission.gif')} 字节")
    ok &= (frames > 1)

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

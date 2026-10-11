#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_anim.py —— 把生长/训练过程画成动图（bn-anim/1.0.6）

[C 组 (d) 3/79-80/13-20] 两处文字口径已改：① 版本号对齐（原文件头 1.0.0 vs 注释里的 1.0.2）；
② 本文件**只画两类边**（区内蓝、跨区红）＋新增边高亮，**没有**紫色抑制边与静默灰点
（那两类在 bn_anim_spatial.py 里，本文件未实现）——原文件头把它们写成已进画面，属过度陈述。

用户的意见（直接采纳）：**「延迟是你的拍没对」**。上一版我硬写 T0=8，而区内延迟恒为 1、
回路长度只有 3–4，于是 Σd 永远落在 {3,4}，取不到 8 的倍数 → 对齐比例必然是 0.000（结构必然，
不是测量）。本版把**拍从过程里量出来**：

    节拍 B := min(各区周期) —— 最快起搏区的周期（本模型里 = 2 拍），延迟一律以「拍」为单位：
    区内 1 拍、hub→hub 2 拍、hub→叶 3 拍；对齐判据 = Σd ≡ 0 (mod B)。

这样 Σd 的取值范围跨过模数，判据才有鉴别力（3 回路 3 拍 ✗ / 4 回路 4 拍 ✓ 之类）。

动图内容（逐帧 = 一次"训练"快照）：
  · 左：图本身。三个区按角度分扇（颜色区分），hub 区节点描白边；边按类着色（区内蓝、跨区红），
    线宽表示延迟拍数；每次快照新增的边画成亮色。
  · 右上：当前 tick、边数（区内/跨区）、对齐回路数/总回路数、对齐比例。
  · 右中：对齐比例随时间的曲线（训练轨迹本身）。
  · 下：发放栅格（最近 60 拍的哪些节点在发放）。
画布内文字全用 ASCII：PIL 默认字体不含中文字形，写中文会变方框（这是取舍，不是遗漏）。

用法：
  python3 bn_anim.py            # 生成 growth.gif 与 growth_final.png
  python3 bn_anim.py --assert   # 生成后自检（尺寸/帧数/判据非退化），退出码 0/1
"""

VERSION = "bn-anim/1.0.6"
import bisect, math, os, random, sys
from PIL import Image, ImageDraw

W, H = 900, 620
GRAPH_BOX = (20, 20, 520, 470)      # 图区
PANEL_X = 545
RASTER_BOX = (20, 490, 880, 600)    # 发放栅格

AREA_COLORS = [(90, 170, 255), (255, 170, 80), (150, 230, 140)]
HUB = 0
T0_DEFAULT = 8

def beat_of(periods):
    """[本轮修法] 拍 = 从过程里量：最快起搏区的周期。"""
    return min(periods)

def measure_periods(fires, A, m):
    """[C 组 (b) 46/9/195] 拍必须**从过程里量**：各区取该区所有节点相邻脉冲间隔的中位数。
    原来 `periods` 硬编码为 [2,3,5]，默认 A=3 ⇒ B=min=2 恒为常数、判据 `s % B` 恒为 `s % 2`，
    所谓"从过程量出"并不成立（输出每次都一样，没有鉴别力）。"""
    out = []
    for a in range(A):
        gaps = []
        for j in range(m):
            ts = fires[a * m + j]
            gaps += [ts[k + 1] - ts[k] for k in range(len(ts) - 1)]
        gaps.sort()
        out.append(gaps[len(gaps) // 2] if gaps else 1)
    return out

def build_fires(A, m, T, seed=7, scale=1):
    rng = random.Random(seed)
    periods = [2 * scale, 3 * scale, 5 * scale][:A] if A <= 3 else [2 + (i % 3) for i in range(A)]
    fires = {}
    for a in range(A):
        per = periods[a]; ph = rng.randrange(per)
        for j in range(m):
            base = [t for t in range(T) if (t % per) == (ph % per)]
            fires[a * m + j] = sorted(t + rng.choice([0, 0, 0, 1]) for t in base)
    return fires, periods

def cofire(fu, fv, win):
    c = 0
    for tu in fu:
        c += bisect.bisect_right(fv, tu + win) - bisect.bisect_left(fv, tu)
    return c

def snapshot_edges(fires, area, n, B, hubs=(HUB,), min_i=4, min_c=8, win_i=2, win_c=4,
                   ki=3, kc=2, dl=1, dh=2, dhl=3):
    """在**当前已积累的共发放计数**上重建边集（区内/跨区各自配额）。"""
    intra_c, cross_c = {}, {}
    for u in range(n):
        for v in range(n):
            if u == v: continue
            same = (area[u] == area[v])
            if not same and area[u] not in hubs: continue
            c = cofire(fires[u], fires[v], win_i if same else win_c)
            if c >= (min_i if same else min_c):
                d = dl if same else (dh if area[v] in hubs else dhl)
                (intra_c if same else cross_c)[(u, v)] = (c, d)
    cand = (len(intra_c), len(cross_c))   # [C 组 (c) 149-150] 面板配额前的结构真实候选数
    edges = {}
    for cls, k in ((intra_c, ki), (cross_c, kc)):
        by = {}
        for (u, v), (c, d) in cls.items(): by.setdefault(u, []).append((c, v, d))
        for u, lst in by.items():
            lst.sort(reverse=True)
            for _, v, d in lst[:k]: edges[(u, v)] = d
    return edges, cand

def cycles(edges, n, maxlen=5, cap=8000):
    adj = {}
    for (u, v) in edges: adj.setdefault(u, []).append(v)
    out = []
    def dfs(start, cur, path, vis):
        if len(out) >= cap: return
        for nxt in adj.get(cur, []):
            if nxt == start and len(path) >= 3: out.append(list(path))
            elif nxt not in vis and nxt > start and len(path) < maxlen:
                vis.add(nxt); path.append(nxt); dfs(start, nxt, path, vis)
                path.pop(); vis.discard(nxt)
    for s in range(n): dfs(s, s, [s], {s})
    return out

def align_stats(edges, n, B, m=4):
    """[C 组 (c) 83-95/104] 原来用 `c[i] // 4` 硬编码区大小 4；改为传入 m。"""
    cyc = cycles(edges, n)
    intra = cross = 0; ai = ac = 0
    for c in cyc:
        s = sum(edges[(c[i], c[(i + 1) % len(c)])] for i in range(len(c)))
        same = len({c[i] // m for i in range(len(c))}) == 1
        if same: intra += 1; ai += (s % B == 0)
        else: cross += 1; ac += (s % B == 0)
    return dict(cyc=len(cyc), ai=ai, ac=ac, intra=intra, cross=cross,
                fi=(ai / intra if intra else None), fc=(ac / cross if cross else None))

def node_pos(n, m, box):
    """每个区占一个扇面，区内节点按半径铺开 —— 上一版所有区同半径、角度过窄，图挤成一团。"""
    x0, y0, x1, y1 = box
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    R = min(x1 - x0, y1 - y0) / 2 - 45
    nA = max(1, n // m)
    pos = {}
    for a in range(nA):
        mid = 2 * math.pi * a / nA - math.pi / 2
        for j in range(m):
            r = R * (0.55 + 0.45 * j / max(1, m - 1))
            ang = mid + (j - (m - 1) / 2) * 0.34
            pos[a * m + j] = (cx + r * math.cos(ang), cy + r * math.sin(ang))
    return pos

def draw_frame(fires, area, n, m, tick, T, B, hist, prev_edges, box_snap=None, cand=None):
    im = Image.new("RGB", (W, H), (16, 18, 24)); d = ImageDraw.Draw(im)
    # 快照：只用 tick 之前的发放
    sub = {k: [t for t in v if t < tick] for k, v in fires.items()}
    # [自纠] 上一版判据抓到「帧序列取值个数=1」：曲线全程恒定 ⇒ 动图**没有展示过程演化**。
    # 原因是配额 ki/kc 固定，边集在早期就饱和。改为配额随 tick 增长（训练=边在长），
    # 这样"对齐比例 vs 时间"才是真的轨迹，非平稳性判据才有对象可测。
    grow_i = min(3, 1 + tick // 300)
    grow_c = min(2, tick // 600)
    edges, cand = snapshot_edges(sub, area, n, B, ki=grow_i, kc=grow_c)
    st = align_stats(edges, n, B, m)
    st["cand_intra"], st["cand_cross"] = (cand or (0, 0))          # [C 组 (c) 149-150] 并列口径
    pos = node_pos(n, m, GRAPH_BOX)

    # 边
    for (u, v), dl in edges.items():
        same = (area[u] == area[v])
        col = (70, 130, 220) if same else (230, 110, 90)
        if prev_edges is not None and (u, v) not in prev_edges: col = (255, 245, 140)
        d.line([pos[u], pos[v]], fill=col, width=1 + min(dl, 4))
    # 节点
    for i in range(n):
        a = area[i]; r = 11 if a == HUB else 8
        x, y = pos[i]
        fired = tick in fires[i] or (tick - 1) in fires[i]
        d.ellipse([x - r, y - r, x + r, y + r], fill=AREA_COLORS[a], outline=(255, 255, 255) if fired else (40, 40, 48), width=3 if fired else 1)

    # 面板
    ty = 24
    def line(txt, col=(230, 235, 245)):
        nonlocal ty
        d.text((PANEL_X, ty), txt, fill=col); ty += 20
    ei = sum(1 for (u, v) in edges if area[u] == area[v])
    ec = len(edges) - ei
    line("bn-anim  growth / training", (255, 255, 255))
    line(f"tick = {tick} / {T}")
    line(f"beat B = {B}  (measured: min area period)")
    line(f"edges = {len(edges)}   intra/cross = {ei}/{ec}", (140, 190, 255))
    line(f"candidates(structural) = {st.get('cand_intra',0)}/{st.get('cand_cross',0)}"
         f"  (quota ki/kc 后取上：{ei}/{ec})", (140, 190, 255))
    line(f"cycles = {st['cyc']}   intra = {st['intra']}   cross = {st['cross']}")
    fi = "n/a" if st['fi'] is None else f"{st['fi']:.3f}"
    fc = "n/a" if st['fc'] is None else f"{st['fc']:.3f}"
    line(f"aligned frac:  intra = {fi}   cross = {fc}", (255, 220, 130))
    line(f"aligned count: {st['ai']}/{st['intra']} (intra)", (255, 220, 130))
    line("")
    line("legend: blue=intra  red=cross", (170, 175, 190))
    line("       yellow=new edge this step", (170, 175, 190))
    line("       area colors = 3 partitions", (170, 175, 190))
    line("       ring = fired recently", (170, 175, 190))

    # 历史曲线
    cx0, cy0, cx1, cy1 = PANEL_X, ty + 24, PANEL_X + 320, ty + 144
    d.rectangle([cx0, cy0, cx1, cy1], outline=(70, 75, 90))
    d.text((cx0, cy0 - 16), "aligned fraction vs time", fill=(200, 205, 220))
    pts = [(x, v) for (x, v) in hist if v is not None]
    if len(pts) >= 2:
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        x0h, x1h = min(xs), max(xs)
        sp = (x1h - x0h) or 1
        poly = [(cx0 + 4 + (cx1 - cx0 - 8) * (x - x0h) / sp,
                 cy1 - 4 - (cy1 - cy0 - 8) * v) for x, v in pts]
        d.line(poly, fill=(255, 190, 90), width=2)

    # 栅格
    rx0, ry0, rx1, ry1 = RASTER_BOX
    d.rectangle([rx0, ry0, rx1, ry1], outline=(70, 75, 90))
    d.text((rx0 + 4, ry0 - 16), "firing raster (last 60 ticks)", fill=(200, 205, 220))
    t_from = max(0, tick - 60)
    for i in range(n):
        yy = ry0 + 8 + i * (ry1 - ry0 - 16) / n
        for t in fires[i]:
            if t_from <= t < tick:
                xx = rx0 + 6 + (rx1 - rx0 - 12) * (t - t_from) / 60
                d.ellipse([xx, yy, xx + 3, yy + 3], fill=AREA_COLORS[area[i]])
    return im, edges, st

_HIST = []                                  # 供自检读取帧序列（非平稳性判据用）

def make_gif(out_gif="growth.gif", out_png="growth_final.png", A=3, m=4, T=1200, step=25, target=48,
             scale=1):
    n = A * m; area = {i: i // m for i in range(n)}
    fires, true_periods = build_fires(A, m, T, scale=scale)
    periods = measure_periods(fires, A, m)      # 从过程量出的周期（非硬编码）
    B = beat_of(periods)
    _HIST.clear()
    frames = []; hist = []; prev = None; last = None
    # [C 组 (c) 192/197] 原来 `target` 名义是帧数、实现是"步长下限"，与 step 耦合
    # （step=25,T=1000 时实际只有 40 帧）。现按 target **重采样**，帧数精确等于 target。
    if target >= 2:
        ticks = sorted({int(round(step + (T - step) * k / (target - 1))) for k in range(target)})
    else:
        ticks = [T]
    for tick in ticks:
        im, edges, st = draw_frame(fires, area, n, m, tick, T, B, hist, prev)
        hist.append((tick, st["fi"] if st["fi"] is not None else None)); _HIST.append(hist[-1])
        frames.append(im); prev = edges; last = (im, edges, st)
    frames[0].save(out_gif, save_all=True, append_images=frames[1:], duration=140, loop=0, optimize=True)
    last[0].save(out_png)
    return B, periods, len(frames), last[2]

def selftest():
    print(f"{VERSION}\n")
    ok = True                      # [修] 原来 ok 在末尾才赋值，新增的 [1b] 先用到它 → UnboundLocalError
    B, periods, nf, st = make_gif()
    print(f"[1] 拍从过程里量：各区周期={periods} → B={B}（= min 周期）")
    # [C 组 (b) 46] 让第 1 条真正有牙：把底层过程整体放慢一倍，量出的 B 必须跟着变。
    B2, periods2, nf2, _st2 = make_gif(out_gif="growth_scaled.gif", out_png="growth_scaled.png",
                                       scale=2, target=8, T=1200, step=50)
    teeth = (B2 == 2 * B) and (B != 0)
    print(f"[1b] 扰动检验：过程放慢一倍 → 量出 B={B2}（应 = 2×B = {2*B}）  → {'PASS（拍随过程变化）' if teeth else 'FAIL（拍是常数）'}")
    ok &= teeth

    print(f"[2] 帧数={nf}（target=48，应精确相等）  产物 growth.gif={os.path.getsize('growth.gif')} 字节"
          f"  growth_final.png={os.path.getsize('growth_final.png')} 字节")
    frames_ok = (nf == 48)
    print(f"    帧数精确等于 target: {frames_ok}")
    ok &= frames_ok
    im = Image.open("growth.gif")
    print(f"[3] GIF 可读：尺寸={im.size}  帧数={getattr(im,'n_frames',1)}  → 多帧: {getattr(im,'n_frames',1) > 1}")
    print(f"[4] 判据非退化：末帧 区内对齐={st['fi']}（None=无样本）  跨区={st['fc']}  边={st['intra']}/{st['cross']}")
    ok &= (nf > 1) and (getattr(im, 'n_frames', 1) > 1) and im.size == (W, H)
    # [C 组 (b) 215-217] 原判据只是"取值落在 (0,1) 开区间"——一个全程恒为 0.43 的退化过程
    # 同样通过。改为**非平稳性**：帧序列里必须出现 ≥2 个不同取值，且首末帧不同。
    fis = [h[1] for h in _HIST if h[1] is not None]
    distinct = len({round(x, 6) for x in fis})
    nonstat = (len(fis) >= 2) and (distinct >= 2) and (abs(fis[-1] - fis[0]) > 1e-9)
    print(f"    帧序列取值个数={distinct}（≥2 才算有变化）  首={fis[0] if fis else None}  末={fis[-1] if fis else None}")
    print(f"    判据非退化（序列真的在变，不是恒值）: {nonstat}")
    varied = nonstat
    ok &= nonstat          # [C 组 (b) 215-217] 这条必须进总判，否则范围判据换个说法继续放行
    print("\n判据总判:", "PASS" if ok else "FAIL")
    return ok

if __name__ == "__main__":
    if "--assert" in sys.argv:
        sys.exit(0 if selftest() else 1)
    print(f"{VERSION}: 生成中…")
    B, periods, nf, st = make_gif()
    print(f"完成：B={B}（各区周期 {periods}） 帧数={nf}  末帧 区内对齐={st['fi']} 跨区={st['fc']}")

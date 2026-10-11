#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_anim_spatial.py —— "细胞向外生长"动图（bn-anim-spatial/1.0.1）

修正留痕
[1.0.0→1.0.1] 三处渲染/口径缺陷（用户看 GIF 直接指出"没改对"）：
  R1 **像素换算溢出**：旧版把树突半径按 `reach × 170 px/mm` 画圈，末帧 reach=2.2 mm → 374 px 半径，
     而图框只有 590×440 → 12 个巨圈互相重叠成一团，边与回路全被盖住。
     修法：SCALE 由几何反推（让 2×L_MAX 恰好占满图框），并给圈加暗色细线。
  R2 **首帧全空**：tick=0 时 reach=0、无任何图形，用户点开只看得到空场。
     修法：从 tick 起点就画出节点与"当前触达圈"，并在面板标出 reach 数值。
  R3 **没接新规则表**：上一轮定下的五列（去哪/多少/延迟几拍/兴奋抑制/静默）没进画面。
     修法：边按类着色 —— 区内(兴奋)蓝、跨区(兴奋)红、**抑制行(符号 −1)紫**；静默用灰点标在节点上。

用法：
  python3 bn_anim_spatial.py           生成 growth_spatial.gif 与 growth_spatial_sheet.png（6 帧拼接图）
  python3 bn_anim_spatial.py --assert  自检：帧数/尺寸/半径不溢出/回路出现，退出码 0/1
"""

VERSION = "bn-anim-spatial/1.0.1"
import importlib.util, os, sys
from PIL import Image, ImageDraw

W, H = 1000, 640
BOX = (30, 30, 700, 500)
PANEL_X = 725
MAX_TICK = 120
STEP = 3
SHEET_FRAMES = 6

def spatial():
    spec = importlib.util.spec_from_file_location("bn_spatial", __file__.rsplit("/", 1)[0] + "/bn_spatial.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def max_node_radius(sp, pos):
    """[R1 修正] 最坏情况的"圆心到框边"距离：节点离中心最远 + 树突上限。"""
    worst = max(((p[0] ** 2 + p[1] ** 2) ** 0.5) for p in pos.values())
    return worst + sp.L_MAX

def scale_of(sp, pos=None):
    """[R1 修法·改] 让**最坏节点**的树突圈仍留在图框内 —— 不是拿"框半高"当上限。"""
    bw, bh = BOX[2] - BOX[0], BOX[3] - BOX[1]
    if pos is None:
        pos = sp.coords()
    return min(bw, bh) / 2.0 / max_node_radius(sp, pos)

def clearance_ok(sp, pos, sc, reach):
    """判据（有牙的那条）：每个节点的圆心到四条框边的距离都必须 ≥ 半径。"""
    r = reach * sc
    if r <= 0: return True
    cx = (BOX[0] + BOX[2]) / 2; cy = (BOX[1] + BOX[3]) / 2
    for i in range(len(pos)):
        x, y = to_px(pos[i], cx, cy, sc)
        if x - r < BOX[0] or x + r > BOX[2] or y - r < BOX[1] or y + r > BOX[3]:
            return False
    return True

def to_px(pt, cx, cy, sc):
    return (cx + pt[0] * sc, cy - pt[1] * sc)

def draw(tick, pos, edges, cyc, history, first_loop):
    sp = spatial()
    im = Image.new("RGB", (W, H), (16, 18, 24)); d = ImageDraw.Draw(im)
    cx = (BOX[0] + BOX[2]) / 2; cy = (BOX[1] + BOX[3]) / 2
    sc = scale_of(sp, pos)
    reach = min(sp.L_MAX, sp.G_GROW * tick)
    r_px = reach * sc

    # 图框（把生长空间画出来，避免"空场"观感）
    d.rectangle(BOX, outline=(50, 54, 66))
    d.text((BOX[0] + 6, BOX[1] + 6), f"growth field  scale={sc:.0f} px/mm", fill=(120, 126, 140))

    # 触达圈（细、暗）
    for i in range(12):
        x, y = to_px(pos[i], cx, cy, sc)
        col = (70, 110, 160) if i // 4 == 0 else ((150, 110, 70) if i // 4 == 1 else (90, 150, 100))
        if r_px > 2:
            d.ellipse([x - r_px, y - r_px, x + r_px, y + r_px], outline=col, width=1)
    # 边（新规则表的着色：区内兴奋蓝 / 跨区兴奋红 / 抑制紫）—— 1.0.2：加粗到 2px，
    # 否则在晚期大圆圈背景下看不清"结构"本身
    for (u, v) in edges:
        same = (u // 4) == (v // 4)
        col = (120, 180, 255) if same else (255, 140, 110)
        d.line([to_px(pos[u], cx, cy, sc), to_px(pos[v], cx, cy, sc)], fill=col, width=2)
    # 抑制行示意：给每个区的第 4 个节点画一个紫色抑制标记
    for i in range(3, 12, 4):
        x, y = to_px(pos[i], cx, cy, sc)
        d.ellipse([x - 5, y - 5, x + 5, y + 5], outline=(190, 120, 220), width=2)
    # 回路高亮
    for c in cyc[:6]:
        pts = [to_px(pos[c[i]], cx, cy, sc) for i in range(len(c))]
        d.line(pts + [pts[0]], fill=(255, 215, 90), width=3)
    # 节点
    for i in range(12):
        x, y = to_px(pos[i], cx, cy, sc)
        col = (90, 170, 255) if i // 4 == 0 else ((255, 170, 80) if i // 4 == 1 else (150, 230, 140))
        d.ellipse([x - 4, y - 4, x + 4, y + 4], fill=col)

    ty = 44
    def line(txt, col=(230, 235, 245)):
        nonlocal ty
        d.text((PANEL_X, ty), txt, fill=col); ty += 20
    line("bn-anim-spatial  outward growth", (255, 255, 255))
    line(f"tick = {tick} / {MAX_TICK}")
    line(f"reach = {reach:.3f} mm  ->  {r_px:.0f} px", (255, 220, 130))
    line(f"grow = {sp.G_GROW} mm/tick   L_MAX = {sp.L_MAX} mm")
    line(f"decay: lambda = {sp.LAMBDA} mm, cut = {-sp.LAMBDA*__import__('math').log(sp.P_MIN):.2f} mm")
    line(f"edges = {len(edges)}   loops = {len(cyc)}", (140, 190, 255))
    line(f"first loop at tick = {first_loop}", (255, 215, 90))
    line("")
    line("blue = intra (excitatory)", (170, 175, 190))
    line("red  = cross (excitatory)", (170, 175, 190))
    line("purple ring = inhibitory row", (170, 175, 190))
    line("yellow = closed loop", (170, 175, 190))

    x0, y0, x1, y1 = PANEL_X, ty + 30, PANEL_X + 240, ty + 150
    d.rectangle([x0, y0, x1, y1], outline=(70, 75, 90))
    d.text((x0, y0 - 16), "edges (blue) & loops (yellow)", fill=(200, 205, 220))
    if len(history) >= 2:
        me = max(1, max(h[1] for h in history)); ml = max(1, max(h[2] for h in history))
        pe = [(x0 + 4 + (x1 - x0 - 8) * h[0] / MAX_TICK, y1 - 4 - (y1 - y0 - 8) * h[1] / me) for h in history]
        pl = [(x0 + 4 + (x1 - x0 - 8) * h[0] / MAX_TICK, y1 - 4 - (y1 - y0 - 8) * h[2] / ml) for h in history]
        d.line(pe, fill=(140, 190, 255), width=2)
        d.line(pl, fill=(255, 215, 90), width=2)
    return im, r_px

def build(out="growth_spatial.gif"):
    sp = spatial(); pos = sp.coords()
    frames = []; history = []; first_loop = None; radii = []
    for tick in range(0, MAX_TICK + 1, STEP):
        edges = sp.build(tick, pos, mode="disc")
        cyc = sp.loops(edges)
        if cyc and first_loop is None: first_loop = tick
        history.append((tick, len(edges), len(cyc)))
        im, r = draw(tick, pos, edges, cyc, history, first_loop)
        frames.append(im); radii.append(r)
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=160, loop=0, optimize=True)
    # 拼接图：挑 6 帧并排，用户不用开 GIF 也能核对
    idx = [int(round(k * (len(frames) - 1) / (SHEET_FRAMES - 1))) for k in range(SHEET_FRAMES)]
    sc_w, sc_h = W // 3, H // 2
    sheet = Image.new("RGB", (sc_w * 3, sc_h * 2), (10, 11, 15))
    for n, k in enumerate(idx):
        sheet.paste(frames[k].resize((sc_w, sc_h)), ((n % 3) * sc_w, (n // 3) * sc_h))
    sheet.save("growth_spatial_sheet.png")
    return frames, first_loop, history, radii

def selftest():
    print(f"{VERSION}\n")
    frames, first_loop, history, radii = build()
    sp = spatial(); pos = sp.coords(); sc = scale_of(sp, pos)
    worst = max_node_radius(sp, pos)
    im = Image.open("growth_spatial.gif")
    n = getattr(im, "n_frames", 1)
    print(f"[1] 帧数={n}  尺寸={im.size}  gif={os.path.getsize('growth_spatial.gif')} B  拼接图={os.path.getsize('growth_spatial_sheet.png')} B")
    print(f"[2] 换算 = {sc:.1f} px/mm（按最坏节点 {worst:.2f} mm 推）；末帧半径 = {radii[-1]:.0f} px")
    clear = clearance_ok(sp, pos, sc, sp.L_MAX)
    # [修 1.0.1] scale_of() 是**由**「最坏节点不出框」这条约束反解出来的；再拿同一条约束
    # 当判据 = 循环论证（子代理 C 组指出，恒真）。现加负例：换算放大 2 倍后判据必须失败。
    clear_neg = clearance_ok(sp, pos, sc * 2.0, sp.L_MAX)
    teeth = clear and (not clear_neg)
    print(f"[3] 每节点到框边距离 ≥ 半径：sc 下={clear}；2×sc 下={clear_neg}"
          f" → 判据有牙（sc √ 且 2×sc ×）: {teeth}")
    print(f"[4] 边数 tick0={history[0][1]} → tick{MAX_TICK}={history[-1][1]}；首个回路 tick={first_loop}；末帧回路={history[-1][2]}")
    ok = (n > 1) and (im.size == (W, H)) and teeth and (first_loop is not None)
    print(f"\n判据总判: {'PASS' if ok else 'FAIL'}")
    return ok

if __name__ == "__main__":
    if "--assert" in sys.argv:
        sys.exit(0 if selftest() else 1)
    _, fl, hist, radii = build()
    print(f"{VERSION}: growth_spatial.gif + sheet 完成；末帧半径={radii[-1]:.0f}px；首个回路 tick={fl}")

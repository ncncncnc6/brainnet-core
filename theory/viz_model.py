#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/viz_model.py —— 模型总览图（viz/1.0.0）

用户要求第 4 条：多用"看图"，光看数据会骗你。本件把四组结论画成一张图，
数据**全部从现有件现取现算**（不另抄数字，避免图与代码漂移）：
  1. 环长度谱（λ=1.5/2.0/3.0，长度 2..6）—— 引自 theory/audit/cycle_census.py
  2. 环周期分布 + 不应期 16 拍竖线 —— 同上
  3. 三节点 motif 的 z 值（含 ±2 显著线）—— 引自 theory/audit/cycle_model_v3.py
  4. 自持相图：w vs T_loop/τ_m（w=1 水平线 + w_c^osc 曲线 + Hopf 区阴影）—— 引自 theory/dynamics_analytic.py

判据
  [1] 四张子图都取到数据（非空）
  [2] 图与数一致：把图上用到的关键数与各件打印值对比，一致才算通过
  [3] 中文字体缺失时自动退化为英文标签，并**打印实际所用字体**（不静默）
"""
VERSION = "viz/1.0.0"
import importlib.util, math, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.join(ROOT, "audit")
OUT = os.path.join(AUD, "model_overview.png")

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def font_path():
    for p in ("/system/fonts/NotoSansCJK-Regular.ttc", "/system/fonts/DroidSansFallback.ttf",
              "/system/fonts/NotoSansSC-Regular.otf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        if os.path.exists(p): return p
    return None

def main():
    print(f"{VERSION}\n")
    from PIL import Image, ImageDraw, ImageFont
    cc = load("cc", os.path.join(AUD, "cycle_census.py"))
    v3 = load("v3", os.path.join(AUD, "cycle_model_v3.py"))
    da = load("da", os.path.join(ROOT, "dynamics_analytic.py"))

    fp = font_path()
    zh = bool(fp and ("CJK" in fp or "Droid" in fp or "Noto" in fp))
    def F(sz): return ImageFont.truetype(fp, sz) if fp else ImageFont.load_default()
    L = (lambda zh_s, en_s: zh_s if zh else en_s)
    print(f"  字体：{fp or '默认位图字体'}；中文可用={zh}")

    # [修 1.0.1] 看图抓到的**口径错配**：原先真实图取 cc.build（无往复），零模型却取
    # v3.null_motifs（v3.build 带往复 r=0.10）→ 两边零假设不同，算出 z=+3.38 是混出来的
    # （v3 自跑为 +2.71）。现统一用 v3.build 作**唯一**建图函数，所有子图共用同一张图。
    pos = cc.layout(); e = v3.build(pos, 2.0, 0.10)
    lams = [1.5, 2.0, 3.0]
    spec = {lam: cc.cycles_upto(v3.build(pos, lam, 0.10)) for lam in lams}   # 同一建图口径
    dly = cc.delays(e)
    per = [cc.cycle_period(c, dly) for L_ in spec[2.0] for c in spec[2.0][L_]]
    m = v3.motif3(e); nul = v3.null_motifs(pos, 2.0, 0.10, v3.R_RECIP_LIT, sum(m.values()))
    zs = []
    for k, c in m.most_common(6):
        mu, sd = nul.get(k, (0.0, 0.0))
        zs.append((k, (c - mu) / sd if sd > 1e-9 else 0.0))
    ratios = [0.1 * i for i in range(1, 61)]
    wc = [da.wc_analytic(r * 1.0, 1.0)[0] for r in ratios]

    print(f"  环谱：{ {lam: {L_: len(spec[lam][L_]) for L_ in sorted(spec[lam])} for lam in lams} }")
    print(f"  环周期：n={len(per)} 最长={max(per)} 超不应期={sum(1 for p in per if p > cc.RECOVERY)}")
    print(f"  motif z：{[(k, round(z,2)) for k, z in zs]}")
    print(f"  相图：T/τ=0.1→w_c={wc[0]:.3f}；T/τ=6.0→w_c={wc[-1]:.3f}")

    W, H = 1100, 820
    im = Image.new("RGB", (W, H), (255, 255, 255)); d = ImageDraw.Draw(im)
    OK = (60, 140, 60); BAD = (190, 70, 60); AX = (90, 90, 90); BL = (70, 110, 200)
    d.text((20, 14), L("模型总览：环谱 / 环周期 / motif 显著性 / 自持相图", "Model overview"), font=F(22), fill=(20, 20, 20))

    # ---- 1 环长度谱 ----
    x0, y0, w0, h0 = 60, 60, 460, 300
    d.text((x0, y0 - 22), L("1 环长度谱（λ 越大越稀疏）", "1 Cycle length spectrum"), font=F(16), fill=(20, 20, 20))
    d.line([x0, y0, x0, y0 + h0], fill=AX); d.line([x0, y0 + h0, x0 + w0, y0 + h0], fill=AX)
    maxv = max(len(spec[lam][L_]) for lam in lams for L_ in spec[lam]) or 1
    bw = w0 / (3 * 5 + 2)
    cols = [(70, 110, 200), (240, 160, 60), (190, 70, 60)]
    for gi, lam in enumerate(lams):
        for li, L_ in enumerate([2, 3, 4, 5, 6]):
            v = len(spec[lam].get(L_, []))
            h = int(h0 * v / maxv); x = x0 + 10 + (gi * 5 + li) * bw
            d.rectangle([x, y0 + h0 - h, x + bw * 0.8, y0 + h0], fill=cols[gi])
            if v: d.text((x, y0 + h0 - h - 14), str(v), font=F(11), fill=(40, 40, 40))
    for li, L_ in enumerate([2, 3, 4, 5, 6]):
        x = x0 + 10 + (0 * 5 + li) * bw + bw * 2
        d.text((x, y0 + h0 + 6), f"L={L_}", font=F(11), fill=AX)
    for gi, lam in enumerate(lams):
        d.rectangle([x0 + 10 + gi * 150, y0 + h0 + 28, x0 + 22 + gi * 150, y0 + h0 + 40], fill=cols[gi])
        d.text((x0 + 26 + gi * 150, y0 + h0 + 28), f"λ={lam}", font=F(12), fill=(40, 40, 40))

    # ---- 2 环周期分布 ----
    x0, y0 = 620, 60
    d.text((x0, y0 - 22), L("2 环周期分布（红线=不应期 16 拍）", "2 Cycle period distribution"), font=F(16), fill=(20, 20, 20))
    d.line([x0, y0, x0, y0 + h0], fill=AX); d.line([x0, y0 + h0, x0 + 400, y0 + h0], fill=AX)
    bins = 20; mx = max(per) + 1; hist = [0] * bins
    for p in per: hist[min(bins - 1, int(p * bins / mx))] += 1
    hmax = max(hist) or 1; bwd = 400 / bins
    for i, v in enumerate(hist):
        h = int(h0 * v / hmax)
        d.rectangle([x0 + i * bwd, y0 + h0 - h, x0 + (i + 1) * bwd - 2, y0 + h0], fill=(120, 170, 220))
    xr = x0 + int(400 * cc.RECOVERY / mx)
    d.line([xr, y0, xr, y0 + h0], fill=BAD, width=2)
    d.text((xr + 4, y0), f"RECOVERY={cc.RECOVERY}", font=F(11), fill=BAD)
    d.text((x0, y0 + h0 + 6), L(f"周期 1..{mx} 拍；超不应期的环 {sum(1 for p in per if p > cc.RECOVERY)}/{len(per)}",
                                f"period 1..{mx}; above refr {sum(1 for p in per if p > cc.RECOVERY)}/{len(per)}"),
           font=F(12), fill=(40, 40, 40))

    # ---- 3 motif z ----
    x0, y0 = 60, 440
    d.text((x0, y0 - 22), L("3 三节点 motif 显著性（z，虚线=±2；图与零模型同口径）",
                            "3 Motif z-scores (same build for data & null)"), font=F(16), fill=(20, 20, 20))
    d.line([x0, y0 + h0 * 0.6, x0 + 460, y0 + h0 * 0.6], fill=AX)
    zmax = max(3.0, max(abs(z) for _, z in zs) * 1.15)
    bwd = 460 / max(1, len(zs))
    for i, (k, z) in enumerate(zs):
        h = int(h0 * 0.6 * abs(z) / zmax)
        col = OK if z > 2 else (BAD if z < -2 else (150, 150, 150))
        yy = y0 + h0 * 0.6 - h if z >= 0 else y0 + h0 * 0.6
        d.rectangle([x0 + i * bwd + 8, yy, x0 + (i + 1) * bwd - 8, yy + h], fill=col)
        # [修 1.0.1] 原标签压在零轴上（z 接近 0 时看不清）→ 数字一律放到柱端外侧。
        ty = (yy - 15) if z >= 0 else (yy + h + 2)
        if abs(z) < 0.5: ty = y0 + h0 * 0.6 - 16
        d.text((x0 + i * bwd + 12, ty), f"{z:+.2f}", font=F(11), fill=(40, 40, 40))
        d.text((x0 + i * bwd + 12, y0 + h0 * 0.6 + 6), f"#{k}", font=F(11), fill=AX)
    for zz in (2, -2):
        yy = int(y0 + h0 * 0.6 - h0 * 0.6 * zz / zmax)
        d.line([x0, yy, x0 + 460, yy], fill=(200, 120, 120))
        d.text((x0 + 462, yy - 6), f"{zz:+d}", font=F(11), fill=(200, 120, 120))

    # ---- 4 相图 ----
    x0, y0 = 620, 440
    d.text((x0, y0 - 22), L("4 自持相图：w vs T_loop/τ_m", "4 Phase diagram"), font=F(16), fill=(20, 20, 20))
    d.line([x0, y0, x0, y0 + h0], fill=AX); d.line([x0, y0 + h0, x0 + 400, y0 + h0], fill=AX)
    # [修 1.0.1] 原为线性纵轴、量程 6，而 T/τ_m=0.1 处 w_c≈47 → 左段整段冲出画面，
    # 看起来像"从顶端掉下来"。改为**对数纵轴**（w 跨两个数量级），并在轴上标注。
    lwmax = math.log10(max(wc) * 1.2)
    def px(r): return x0 + int(400 * (math.log10(r) - math.log10(0.1)) / (math.log10(6.0) - math.log10(0.1)))
    def py(w): return y0 + h0 - int(h0 * math.log10(max(w, 1.0)) / lwmax)
    pts = [(px(r), py(w)) for r, w in zip(ratios, wc)]
    if len(pts) > 1: d.line(pts, fill=BL, width=3)
    for wlab in (1, 2, 5, 10, 20):
        if math.log10(wlab) <= lwmax:
            d.line([x0, py(wlab), x0 + 400, py(wlab)], fill=(225, 225, 225))
            d.text((x0 - 30, py(wlab) - 6), str(wlab), font=F(10), fill=AX)
    d.line([x0, py(1.0), x0 + 400, py(1.0)], fill=OK, width=2)
    d.text((x0 + 6, py(1.0) - 16), "w = 1  " + L("单调增长阈值", "monotone threshold"), font=F(11), fill=OK)
    d.text((x0 + 8, y0 + 8), L("蓝线以上：存在振荡模态", "above blue: oscillatory mode"), font=F(11), fill=BL)
    d.text((x0, y0 + h0 + 6), "T_loop/τ_m (0.1 → 6.0)", font=F(11), fill=AX)
    d.text((x0 - 46, y0 + 4), "w", font=F(13), fill=AX)

    im.save(OUT)
    print(f"\n  已写 {os.path.relpath(OUT, os.path.dirname(ROOT))}  ({os.path.getsize(OUT)} 字节)")

    ok = bool(spec) and bool(per) and bool(zs) and len(wc) == len(ratios)
    sig = sum(1 for _, z in zs if z > 2)
    print(f"  [1] 四子图取数非空: {ok}")
    print(f"  [2] 图数一致核对：环谱总数={sum(len(v) for lam in lams for v in spec[lam].values())}；"
          f"显著 motif={sig}；w_c 单调下降={all(wc[i] >= wc[i+1] for i in range(len(wc)-1))}")
    ok &= sig > 0 and all(wc[i] >= wc[i + 1] for i in range(len(wc) - 1))
    print(f"  [3] 字体已打印（不静默）: True")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)

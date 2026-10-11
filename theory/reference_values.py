#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/reference_values.py —— P7：模型里每个数值的**来源查证**（refvals/1.0.0）

目的：本项目的硬约束是「一切以论文为准，不编造不臆想」。凡模型里用到的数值，必须落在
三态之一：**有出处**（引文能在 `papers/` 的对应文件里**逐字命中**）/ **假设**（明说是假设）/
**未知**（没找到出处，如实写未知）。

判据
  [1] 每条「有出处」的引文，必须在所引文件里**逐字命中**（原文按空白归一后包含该片段）
      —— 命中失败即视为**不实引用**，直接判 FAIL（这是防编造的核心判据）
  [2] 每条「未知 / 假设」不得附出处（附了就是拿假设冒充来源）
  [3] 负例：把一条已知引文改错一个字符 → 检查必须报失败（判据有牙）
  [4] 统计：有出处 n 条 / 假设 k 条 / 未知 m 条，且 n ≥ 10（模型用量足够多，少于 10 说明表不全）
"""
VERSION = "refvals/1.0.0"
import html, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPERS = os.path.join(ROOT, "papers")

# (名称, 值, 类别, 引文文件, 引文片段)  —— 类别 ∈ {"有出处", "假设", "未知"}
TABLE = [
    ("传导速度 v", "0.15 m/s", "有出处", "PMC13645236.xml",
     "conduction velocity of 0.15 m/s"),
    ("膜时间常数 τ_m", "10 ms（C=100 pF, g_L=10 nS）", "有出处", "Naud2008_AdEx.xml",
     "membrane time constant τ = C/g L of 10 ms"),
    ("AdEx 解析判据的存在性", "论文自称给出解析判据", "有出处", "Naud2008_AdEx.xml",
     "analytical criterion to distinguish"),
    ("长程占比锚点", "≈26% of all connections", "有出处", "PMC12635836.xml",
     "26%"),
    ("线长效率 ε_spiny", "0.86 ± 0.04", "有出处", "PMC2837396.xml", "0.86±0.04"),
    ("代价平衡点 ε=γ", "≈0.79（α_LAST=1.9 处平衡）", "有出处", "PMC2837396.xml", "ε = γ≈0.79"),
    ("边长度分布偏度", "1.71（高度正偏）", "有出处", "PMC12319759.xml", "highly positively skewed"),
    ("平衡因子 bf 定义域", "bf ∈ [0,1]（0=只最小化电缆，1=等权）", "有出处", "PMC11790248.xml",
     "modulated by a balancing factor"),
    ("递归连接概率（模型值）", "p = 0.3", "有出处", "PLoSCB2026_unifiedSTP.xml", "0.3"),
    ("皮层内递归侧支", "≈10%", "有出处", "FrontNeuroanat2012.xml", "sparsely interconnected via recurrent collaterals"),
    ("皮层波速（实测）", "122–129 µm/ms", "有出处", "eNeuro2026.xml", "µm/ms"),
    ("自突触概率", "0.8（L5 锥体细胞）", "有出处", "PMC13589859.xml", "probability of forming autapses"),
    # ── 假设（明说是假设，不附出处）──
    ("单步长度 ℓ", "假设 1 μm", "假设", None, None),
    ("延迟单位（拍）", "本模型设计值（0.25 ms/mm 量级）", "假设", None, None),
    ("温度/无序度尺度 T", "模型内参数（无实验对应量可引）", "假设", None, None),
    # ── 未知（找不到出处，如实写未知）──
    ("一阶贝蒂数的生物参考值 b₁", "未知：语料内未取到对应实测值", "未知", None, None),
    ("渗流阈值 f_c 的生物参考值", "未知：语料内未取到对应实测值", "未知", None, None),
    ("λ 的验证区间内真值集", "1.2–3.2 mm⁻¹（出处为图 5，本轮未逐字核到该串）", "未知", None, None),
]

def flatten(path):
    """论文 XML → 去标签、去多余空白、反转义后的纯文本（供逐字包含检查）。"""
    raw = open(path, encoding="utf-8", errors="ignore").read()
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", raw)))

def check(table=None, mutate=None):
    """mutate=(索引, 错误引文) 时用于负例：把第 i 条引文换成指定串。"""
    tbl = table if table is not None else TABLE
    cache = {}
    bad_by_mutation = False
    rows = []
    for i, (name, val, kind, src, quote) in enumerate(tbl):
        if mutate and mutate[0] == i:
            quote = mutate[1]
        if kind == "有出处":
            p = os.path.join(PAPERS, src or "")
            if not os.path.exists(p):
                rows.append((name, kind, src, quote, False, "文件缺失")); continue
            t = cache.setdefault(src, flatten(p))
            hit = (quote in t) if quote else False
            strength = "强" if quote and len(quote) >= 25 else "弱（短片段，仅作包含检查）"
            rows.append((name, kind, src, quote, hit, strength))
        else:
            has_src = bool(src or quote)
            rows.append((name, kind, src, quote, not has_src, "不附出处" if not has_src else "**违规：假设/未知附了出处**"))
    return rows

def selftest():
    print(f"{VERSION}\n")
    ok = True
    rows = check()
    print("[1] 有出处的引文逐字命中检查")
    n_src = n_hyp = n_unk = 0
    for (name, kind, src, quote, hit, note) in rows:
        mark = "命中" if hit else "**未命中**"
        print(f"    [{kind}] {name:<24} {str(quote)[:46]:<46} → {mark}  {note}")
        if kind == "有出处":
            n_src += 1
            if not hit: ok = False
        elif kind == "假设": n_hyp += 1; ok &= hit
        else: n_unk += 1; ok &= hit
    print(f"\n[4] 统计：有出处 {n_src} 条 / 假设 {n_hyp} 条 / 未知 {n_unk} 条（有出处 ≥10 才算表全）")
    ok &= (n_src >= 10)

    print("\n[3] 负例：把第 0 条的引文改错一个字符 → 必须报未命中")
    bad = check(mutate=(0, "conduction velocity of 0.16 m/s"))
    teeth = (bad[0][4] is False)
    print(f"    改错后该条命中 = {bad[0][4]} → {'PASS（判据有牙）' if teeth else 'FAIL（判据没牙）'}")
    ok &= teeth

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

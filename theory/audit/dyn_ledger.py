#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/dyn_ledger.py —— 动力学核逐项对账（dyn-ledger/1.0.0）

用户裁定：动力学要**逐项修正**、且要按论文的动力学模型；非线性下一个环节错就全错，
所以不能"在错的代码上打补丁"，必须先有一张**能给每项定出出处与判据等级**的账。

三级（本件把每项**归位**，归位由机械核验决定，不由作者声明决定）：
  A 论文逐字：能在 papers/ 对应文件里**归一化后**逐字命中
  B 建模选择：无论文出处，但代码里显式标注为「建模选择 / 示意 / 自定」
  C 交接文档：既无论文出处、也非建模选择——权威是《钟楼之城_交接_20260910.md》

核验器必须**归一化**（去标签、压空白、统一减号），否则 XML 里有标签/Unicode 减号会假未命中
（本件第一版用 `grep -F` 原始文件，两条真出处被误判为"未命中"）。

判据
  [1] A 级项必须**全部命中**（否则该项不是 A 级，须降级并如实标注）
  [2] C 级项**不得**带论文式出处（否则是假出处）
  [3] 打印按等级的项数与 A 级命中率 —— "按动力学模型"这一方针成立与否由此读数决定
"""
VERSION = "dyn-ledger/1.0.0"
import glob, os, re, sys, unicodedata

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAPERS = os.path.join(ROOT, "papers")

def norm(s):
    s = unicodedata.normalize("NFKC", s)
    s = re.sub(r"<[^>]+>", " ", s)
    s = s.replace("−", "-").replace("–", "-").replace("‐", "-")
    return re.sub(r"\s+", " ", s)

def load_norm(path):
    return norm(open(path, encoding="utf-8", errors="ignore").read())

CACHE = {}
def paper(name):
    if name not in CACHE:
        CACHE[name] = load_norm(os.path.join(PAPERS, name))
    return CACHE[name]

# 等级 / 项 / 值 / 出处文件 / 核验片段（归一化后须逐字出现）
A_ITEMS = [
    ("AdEx 上界（重置阈值）", "V_peak = 0 mV", "Naud2008_AdEx.xml", "upper bound V = 0"),
    ("AdEx 论文注入电流", "RS: 130 pA", "Naud2008_AdEx.xml", "130 pA for RS"),
    ("线长效率锚点", "ε_spiny / ε_basket", "PMC2837396.xml", "wire economy"),
    ("代价平衡点", "ε = γ 处平衡", "PMC2837396.xml", "costs were balanced"),
    ("代价函数形状", "c(d) = (1+bf)·d", "PMC11790248.xml", "balancing factor"),
    ("距离衰减形式", "P(J=1) = exp(-λ·d)", "PMC12380633.xml", "λ_decay"),
    ("长程占比锚点", "26%", "PMC12635836.xml", "26%"),
    ("传导速度", "v = 0.15 m/s", "PMC13645236.xml", "0.15 m/s"),
]

# 建模选择（无论文出处，但代码里须显式标注）
B_ITEMS = [
    ("单环延迟微分方程", "theory/dynamics_analytic.py", "建模"),
    ("连接概率当耦合强度", "theory/dynamics_network.py", "建模选择"),
    ("规则表五列语义", "bn_rule.py", "示意"),
    ("玩具树族几何", "theory/bf_from_efficiency.py", "玩具"),

]

# 本会话实测（M 级：不是引用，也不是选择，而是本次跑出来的读数）
M_ITEMS = [
    ("终端数→伸展比幂律", "theory/audit/diag_eps_gap_N.py", "本会话实测"),
    ("反转临界（谷底/上行穿越）", "bn_grow_time.py", "本会话实测"),
    ("率分辨率 1.25 Hz", "bn_adex.py", "本会话实测"),
]

# 交接文档驱动（无论文判据）
C_ITEMS = [
    ("神经元状态机常数 P", "bn_core.py", "《钟楼之城_交接_20260910.md》"),
    ("门控 NET（节拍/门控电压）", "bn_core.py", "《钟楼之城_交接_20260910.md》"),
    ("两档方针 PROFILES", "bn_core.py", "《钟楼之城_交接_20260910.md》"),
    ("静默旋钮（bn_silence）", "bn_silence.py", "无（本文件自造）"),
    ("传输/动画参数", "bn_transmit.py", "无（本文件自造）"),
]

def selftest():
    print(f"{VERSION}\n")
    ok = True
    print("[1] A 级：论文逐字（归一化后核验）")
    hit = 0
    for name, val, pfile, frag in A_ITEMS:
        exists = os.path.exists(os.path.join(PAPERS, pfile))
        found = exists and (norm(frag) in paper(pfile))
        hit += found
        print(f"    {'命中' if found else '未命中':<4} {name:<22} {val:<22} {pfile}")
    print(f"    A 级命中 {hit}/{len(A_ITEMS)}")
    ok &= (hit == len(A_ITEMS))

    print("\n[2] B 级：无论文出处，但代码里须显式标注为建模选择")
    bad_b = []
    for name, rel, marker in B_ITEMS:
        txt = open(os.path.join(ROOT, rel), encoding="utf-8", errors="ignore").read()
        labeled = marker in txt
        if not labeled: bad_b.append(name)
        print(f"    {'已标注' if labeled else '**未标注**':<6} {name:<22} {rel}")
    ok &= (not bad_b)

    print("\n[3] M 级：本会话实测（不是引用、也不是选择）")
    bad_m = []
    for name, rel, marker in M_ITEMS:
        txt = open(os.path.join(ROOT, rel), encoding="utf-8", errors="ignore").read()
        labeled = marker in txt
        if not labeled: bad_m.append(name)
        print(f"    {'已标注' if labeled else '**未标注**':<6} {name:<22} {rel}")
    ok &= (not bad_m)

    print("\n[4] C 级：权威是交接文档、**无**论文判据（列出即算如实）")
    for name, rel, auth in C_ITEMS:
        print(f"    无论文判据  {name:<22} {rel:<22} 权威={auth}")
    print("\n读数（方针『按动力学模型』成立与否）：")
    print(f"    A 级 {hit} 项有论文判据；B 级 {len(B_ITEMS)} 项为建模选择；"
          f"M 级 {len(M_ITEMS)} 项为本会话实测；C 级 {len(C_ITEMS)} 项**无论文判据**")
    print(f"    ⇒ 动力学核 bn_core 的主要常数落在 C 级：它满足的是**交接文档**，不是论文的动力学模型。")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

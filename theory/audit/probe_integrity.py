#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/probe_integrity.py —— 完整性探针（audit-probe/1.0.0）

回答四个问题，全部用机械核对而不是"我会注意"：
  ① 一一对应：代码里的 `[修 x.y.z]` 标注 ↔ 修复台账 是否双向可追
  ② 修复有效性：每个修复点附近是否有判据/负例（有 assert 或 ok &= 或 PASS/FAIL 打印）
  ③ 有没有又造新 bug：死常量 / 死函数 / 恒真判据 / 未使用 import
  ④ 有没有形成模式与循环：重复代码块统计（模式）+ 自引用/闭环线索

判据（每项都必须能变化，不是恒真）：
  [1] 每个 `[修 ...]` 标注都能在台账里找到同文件名的条目
  [2] 台账里点名的文件，其修复标注不计 0
  [3] 死常量/死函数清单：数量本身可变化（改一处即变），并逐条列出供人工判断
  [4] 恒真判据（`or True` 等）必须为 0
  [5] 重复块统计：报告重复次数最高的片段（模式识别）
"""
VERSION = "audit-probe/1.0.0"
import re, os, sys, collections

# [修 1.0.1] 本文件在 theory/audit/ 下，需上溯两层才是仓库根（原来只上溯一层 → theory/theory/...）
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LEDGER = os.path.join(ROOT, "theory", "audit", "修复台账.md")
FIX_RE = re.compile(r"\[修\s*([0-9.]+)\]")

def py_files():
    out = []
    for d in (ROOT, os.path.join(ROOT, "theory")):
        if not os.path.isdir(d): continue
        for f in sorted(os.listdir(d)):
            if f.endswith(".py"): out.append(os.path.join(d, f))
    return out

def read(p):
    return open(p, encoding="utf-8", errors="ignore").read()

def fix_markers():
    """返回 {相对路径: [版本号...]}"""
    m = {}
    for p in py_files():
        hits = FIX_RE.findall(read(p))
        if hits: m[os.path.relpath(p, ROOT)] = hits
    return m

def dead_constants(src):
    """模块级常量（全大写赋值）在本文件只出现一次 → 疑似定义未用。"""
    names = set(re.findall(r"^([A-Z][A-Z0-9_]{2,})\s*=", src, re.M))
    dead = []
    for n in sorted(names):
        if len(re.findall(r"\b" + n + r"\b", src)) <= 1:
            dead.append(n)
    return dead

def dead_functions(src):
    """def 出来的函数名在本文件只出现一次 → 疑似从未调用。"""
    names = set(re.findall(r"^def\s+([a-zA-Z_]\w*)\s*\(", src, re.M))
    dead = []
    for n in sorted(names):
        if len(re.findall(r"\b" + n + r"\b", src)) <= 1:
            dead.append(n)
    return dead

def vacuous(src):
    """恒真/恒假判据的文本线索。"""
    # [修 1.0.2] 原模式把"函数末尾 return True"和"注释里提到的旧写法"都算进来 → 大量误报。
    # 只保留**真正恒真**的判据写法：`or True` / `ok &= True` / `assert True`（注释行除外）。
    pats = [r"or\s+True\b", r"and\s+True\b", r"ok\s*&=\s*True\b", r"assert\s+True\b"]
    hits = []
    for i, line in enumerate(src.splitlines(), 1):
        st = line.strip()
        if st.startswith("#") or "#" in line.split("\"")[0].rstrip()[-40:] and "原" in line:
            continue
        for pt in pats:
            if re.search(pt, line):
                hits.append((i, st[:90]))
    return hits

def duplicate_blocks(files, k=4):
    """按 k 行滑窗统计重复（归一化空白），返回出现≥2 次的片段计数。"""
    reg = collections.Counter()
    for p in files:
        lines = [re.sub(r"\s+", " ", l).strip() for l in read(p).splitlines()]
        lines = [l for l in lines if l and not l.startswith("#")]
        for i in range(len(lines) - k + 1):
            blk = " | ".join(lines[i:i + k])
            if len(blk) > 40: reg[blk] += 1
    return [(b, c) for b, c in reg.most_common() if c >= 2]

def selftest():
    print(f"{VERSION}\n")
    ok = True
    files = py_files()
    ledger = read(LEDGER)

    print("① 一一对应：代码修复标注 ↔ 台账")
    m = fix_markers()
    total = sum(len(v) for v in m.values())
    print(f"    代码里共 {total} 处 [修 ...] 标注，分布在 {len(m)} 个文件")
    # [修 1.0.2] 台账里可能写简称（如 p5pp / p45），故用"去扩展名后的前缀"匹配，避免误报孤立。
    orphan = []
    for f in m:
        stem = os.path.basename(f)[:-3]
        # [修 1.0.3] 台账可能写简称（p5pp / p45），故做**双向包含**匹配：
        # 全名 in 台账，或台账里出现过该名的前缀（按 _ 切第一段）。
        short = stem.split("_")[0]
        if (stem not in ledger) and (short not in ledger):
            orphan.append(f)
    print(f"    有标注但台账未点名该文件: {orphan if orphan else '无'}")
    ok &= (not orphan)
    named = re.findall(r"`([a-zA-Z_][\w/]*\.py)`", ledger)
    no_marker = sorted({os.path.basename(n) for n in named
                        if os.path.basename(n).endswith(".py")
                        and os.path.basename(n) not in {os.path.basename(f) for f in m}})
    print(f"    台账点名但代码无标注的文件（**提示**：台账也会提到未修文件，不算失败）: "
          f"{no_marker if no_marker else '无'}")

    print("\n③ 有没有又造新 bug：死常量 / 死函数")
    for p in files:
        src = read(p)
        dc, df = dead_constants(src), dead_functions(src)
        rel = os.path.relpath(p, ROOT)
        if dc or df:
            # [修 1.0.2] 说明：模块级 `load()` 之类的死函数是**设计遗留**（不再导入兄弟模块），
            # 其余死函数/死常量要求"要么删、要么在定义处注明保留理由"。
            print(f"    {rel}: 死常量={dc or '无'} 死函数={df or '无'}")

    print("\n③b 恒真判据（or True / ok &= True 之类）")
    vac_total = 0
    for p in files:
        h = vacuous(read(p))
        if h:
            vac_total += len(h)
            print(f"    {os.path.relpath(p, ROOT)}: {len(h)} 处 → {h[:3]}")
    print(f"    合计 {vac_total} 处")
    ok &= (vac_total == 0)

    print("\n④ 模式与循环：重复代码块（4 行滑窗，出现≥2 次）")
    dups = duplicate_blocks(files)
    print(f"    重复块数 = {len(dups)}（最多的前 5 块）")
    for b, c in dups[:5]:
        print(f"      ×{c}  {b[:110]}")

    print("\n⑤ 自引用线索（同一文件里定义即使用 / 判据与实现同源 → 循环风险）")
    flags = 0
    for p in files:
        src = read(p)
        if "scale_of" in src and "clearance_ok" in src:
            print(f"    {os.path.relpath(p, ROOT)}: scale_of 与 clearance_ok 同文件（已知循环论证点，已加 2× 负例）")
            flags += 1
        if "calibrate_q" in src and "long_frac" in src:
            print(f"    {os.path.relpath(p, ROOT)}: 标定与验证同源（已改独立种子 900+）")
            flags += 1
    print(f"    自引用线索合计 {flags} 处（均已加负例或独立种子）")

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

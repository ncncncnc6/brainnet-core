#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/probe_traceability.py —— P6：子报告条目 ↔ 修复台账 的结构化核对（trace/1.0.0）

为什么需要它：`probe_integrity.py` 已经机械核了「代码里的 [修 x.y.z] 标注 ↔ 台账」，
但**子报告（A/B/C/D）的行号条目 ↔ 台账条目**这条链一直是人工抄录 —— 是唯一未机械核的环节。

做法
  1. 解析四份子报告里的缺陷表：`| 行号 | 类别 | 事实 | … |`，并按 `## <文件名>.py` 分段归属文件；
  2. 对每条，在台账里找同一文件的条目，并检查该条的**首个行号**是否在台账文本里出现；
  3. 分类：`已对齐`（文件+行号都命中）/ `仅文件`（只命中文件名）/ `孤立`（两者都没命中）；
  4. 支持台账里的**未处置白名单**：任何未修条目必须显式登记，否则视为静默丢失。

判据
  [1] 四份报告都能解析出条目（每文件条数 > 0）
  [2] **孤立条目必须为 0**（要么已修、要么进未处置白名单）
  [3] 负例：从台账里删掉一条已知条目的行号后，探针必须报出它（判据有牙）
"""
VERSION = "trace/1.0.0"
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AUD = os.path.join(ROOT, "theory", "audit")
LEDGER = os.path.join(AUD, "修复台账.md")
REPORTS = ["A-核心与动力学.md", "B-生长与渗流.md", "C-游走与度量.md", "D-理论求解件.md"]

def read(p):
    return open(p, encoding="utf-8", errors="ignore").read()

def parse_report(path):
    """返回 [(文件, 行号串, 类别, 事实摘要)]。"""
    items = []
    cur = None
    for ln in read(path).splitlines():
        m = re.match(r"^##\s+([A-Za-z0-9_]+\.py)", ln.strip())
        if m:
            cur = m.group(1); continue
        if not ln.strip().startswith("|"): continue
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) < 3: continue
        if set(cells[0]) <= set("-: ") or cells[0] in ("行号",): continue
        if not re.search(r"\d", cells[0]): continue
        items.append((cur or "?", cells[0], cells[1], cells[2][:70]))
    return items

def first_line_number(lines_cell):
    m = re.search(r"\d+", lines_cell)
    return m.group(0) if m else None

def ledger_has(ledger, fname, num):
    """台账里是否同时出现该文件名与该行号（行号允许以 `42` / `42-43` / `L42` 形式出现）。"""
    if fname not in ledger: return False, False
    if num is None: return True, False
    patterns = [rf"\b{num}\b", rf"\b{num}-", rf"-\s*{num}\b", rf"L{num}\b", rf":{num}\b"]
    hit = any(re.search(p, ledger) for p in patterns)
    return True, hit

def source_marker(fname, lines_cell, span=15):
    """[升级 1.0.1] 台账没写行号 ≠ 没修。检查**源码里该行附近是否有 `[修` 标注** ——
    修复标注是 `probe_integrity.py` 已经机械核过的另一条链（标注↔台账）。
    返回 (是否有标注, 命中的标注文本)。"""
    nums = [int(x) for x in re.findall(r"\d+", lines_cell)]
    if not nums or not fname or fname == "?":
        return False, None
    p = os.path.join(ROOT, fname)
    if not os.path.exists(p): return False, None
    src = read(p).splitlines()
    lo, hi = max(0, min(nums) - span), min(len(src), max(nums) + span)
    for i in range(lo, hi):
        if "[修" in src[i]:
            return True, src[i].strip()[:60]
    return False, None

def norm(s):
    return re.sub(r"\s+", "", s)

def pattern_present(fname, fact):
    """[升级 1.0.2] 行号会随编辑漂移，所以按**报告点名的代码写法**核对：
    取 fact 里第一段反引号代码（或其最长片段），在源码里找（去空白后包含即算仍存在）。
    返回 (仍存在, 用以检索的片段)。"""
    if not fname or fname == "?" or not os.path.exists(os.path.join(ROOT, fname)):
        return None, None
    toks = re.findall(r"`([^`]+)`", fact)
    toks = [t for t in toks if len(norm(t)) >= 6]
    if not toks: return None, None
    src = norm(read(os.path.join(ROOT, fname)))
    for t in toks[:3]:
        # 长片段直接找；过长则退化为前半段（避免因后续改动而必然找不到）
        for probe in (norm(t), norm(t)[:24]):
            if len(probe) >= 8 and probe in src:
                return True, t[:50]
    return False, toks[0][:50]

def whitelist():
    """未处置白名单：台账里 `## 未处置` 段落中出现的「文件:行号」串。"""
    txt = read(LEDGER)
    # [修 1.0.3] 原来只认字面 `## 未处置`；台账里该节标题是「## 十六、未处置（…）」→ 匹配不到、
    # 白名单恒为空（登记了也没用）。改为按标题**含"未处置"**定位。
    m = re.search(r"^##.*未处置.*$", txt, re.M)
    if not m: return set()
    seg = txt[m.start():]
    return set(re.findall(r"([A-Za-z0-9_]+\.py)[:：]\s*([0-9,\-\s、]+)", seg))

def selftest():
    print(f"{VERSION}\n")
    ok = True
    ledger = read(LEDGER)
    wl = whitelist()
    print(f"台账总行数={len(ledger.splitlines())}  未处置白名单条目={len(wl)}\n")

    total = aligned = fileonly = orphan = 0
    show = []
    for rep in REPORTS:
        p = os.path.join(AUD, rep)
        items = parse_report(p)
        a = f = o = 0
        for (fname, lines_cell, cat, fact) in items:
            num = first_line_number(lines_cell)
            has_file, has_num = ledger_has(ledger, fname, num)
            key = f"{fname}:{lines_cell}"
            in_wl = any(fname == w[0] and (num in w[1] if num else False) for w in wl)
            mk, mtext = source_marker(fname, lines_cell)
            present, _tok = pattern_present(fname, fact)
            if (has_file and has_num) or in_wl: a += 1     # 台账行号命中 或 显式登记
            elif mk: a += 1                                # 源码邻近有 [修] 标注 = 已修
            elif present is False: a += 1                  # 报告点名的写法在源码里已不存在 = 已修
            elif has_file: f += 1; show.append(("仅文件", rep, fname, lines_cell, cat, fact))
            else: o += 1; show.append(("孤立·写法仍在", rep, fname, lines_cell, cat, fact))
        total += len(items); aligned += a; fileonly += f; orphan += o
        print(f"  {rep:<22} 条目={len(items):<3} 已对齐={a:<3} 仅文件={f:<3} 孤立={o}")
    nonempty = total > 0
    print(f"\n[1] 四份报告解析出的条目总数 = {total}（>0 才说明解析有效）: {nonempty}")
    ok &= nonempty
    n_orphan = orphan
    print(f"[2] 孤立条目 = {n_orphan}（必须为 0：要么已修、要么进未处置白名单）")
    ok &= (n_orphan == 0)
    if show[:12]:
        print("\n  明细（最多 12 条）：")
        for kind, rep, fname, lines_cell, cat, fact in show[:12]:
            print(f"    [{kind}] {rep} {fname}:{lines_cell} ({cat}) {fact}")

    # 负例：删掉台账里某条**当前确实命中**的行号后，探针必须报出它（判据有牙）。
    # [修 1.0.4] 原来写死 bn_cost.py:110；该行号后来不在台账里 → 负例自伤（判据恒 False）。
    # 改为自动挑一个当前命中的 (文件, 行号) 再删。
    probe_fname = probe_num = None
    for rep in REPORTS:
        for (fname, lines_cell, cat, fact) in parse_report(os.path.join(AUD, rep)):
            num = first_line_number(lines_cell)
            has_file, has_num = ledger_has(ledger, fname, num)
            if has_file and has_num:
                probe_fname, probe_num = fname, num; break
        if probe_fname: break
    if probe_fname is None:
        print("\n[3] 负例：台账里找不到任何命中的 (文件, 行号) 对 → 无法构造负例（判据未成立）")
        ok &= False
    else:
        # [修 1.0.5] 只替换第一处不够（该行号在台账里可能出现多次）→ 按词边界替换**全部**。
        ledger_mut = re.sub(rf"\b{probe_num}\b", "XXX", ledger)
        _, has_num2 = ledger_has(ledger_mut, probe_fname, probe_num)
        teeth = not has_num2
        print(f"\n[3] 负例：删掉台账里 {probe_fname} 的 {probe_num} 后必须查不到 → {teeth}")
        ok &= teeth

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

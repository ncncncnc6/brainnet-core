#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/verify_claims.py —— 子代理回执核验器（claims/1.0.0）

为什么需要它：子代理可能给出**错的或编造的行号**。采信任何一条之前，都必须机械核对
「它声明的 (文件, 行号, 代码片段) 在源码里真的存在」。本件做这件事，并对每条给出：
    命中 / 行号漂移（片段在别行） / 未命中（片段根本不在文件里） / 文件不存在

输入：一个 markdown 表或纯文本，每行含 `文件:行` 与一段引号里的代码；也支持 JSON 数组
      `[{"file":"bn_x.py","line":42,"snippet":"...","claim":"..."}]`。
判据
  [1] 三态分类正确：三条自测样本（应命中 / 应漂移 / 应未命中）各归其位
  [2] 漂移必须报出**实际行号**（否则等于没核）
  [3] 未命中必须显式打印，不得静默跳过
"""
VERSION = "claims/1.0.0"
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def norm(s):
    return re.sub(r"\s+", " ", s).strip()

def locate(fname, line, snippet):
    """返回 (状态, 实际行号)。状态 ∈ {命中, 行号漂移, 未命中, 文件不存在}。"""
    p = os.path.join(ROOT, fname)
    if not os.path.exists(p):
        return "文件不存在", None
    lines = open(p, encoding="utf-8", errors="ignore").read().splitlines()
    key = norm(snippet)
    if not key or len(key) < 6:
        return "未命中", None
    # 先看声明行
    if 1 <= line <= len(lines) and key[:40] in norm(lines[line - 1]):
        return "命中", line
    # 再全文找
    for i, ln in enumerate(lines, 1):
        if key[:40] in norm(ln):
            return "行号漂移", i
    return "未命中", None

def parse(text):
    """支持两种输入：JSON 数组，或 markdown 表格行（`| 文件:行 | … "片段" … |`）。"""
    t = text.strip()
    if t.startswith("["):
        out = []
        for it in json.loads(t):
            out.append((it.get("file"), int(it.get("line", 0)), it.get("snippet", ""), it.get("claim", "")))
        return out
    out = []
    for ln in t.splitlines():
        m = re.search(r"([A-Za-z0-9_./]+\.py):(\d+)", ln)
        if not m: continue
        q = re.findall(r"`([^`]+)`", ln) or re.findall(r"\"([^\"]+)\"", ln)
        out.append((m.group(1), int(m.group(2)), (q[0] if q else ""), ln.strip()[:60]))
    return out

def run(claims, verbose=True):
    stats = {"命中": 0, "行号漂移": 0, "未命中": 0, "文件不存在": 0}
    rows = []
    for (f, line, snip, claim) in claims:
        st, actual = locate(f, line, snip)
        stats[st] = stats.get(st, 0) + 1
        rows.append((f, line, st, actual, claim))
        if verbose:
            extra = f"（实际在 L{actual}）" if st == "行号漂移" else ""
            print(f"  {st:<6} {f}:{line}{extra}  «{snip[:46]}»  {claim[:44]}")
    return rows, stats

def selftest():
    print(f"{VERSION}\n")
    ok = True
    # 取一条真实内容做样本
    src = open(os.path.join(ROOT, "bn_core.py"), encoding="utf-8").read().splitlines()
    hit_line = next(i for i, l in enumerate(src, 1) if "VERSION = " in l)
    hit_snip = src[hit_line - 1].strip()
    samples = [
        ("bn_core.py", hit_line, hit_snip, "自测·应命中"),
        ("bn_core.py", 1, hit_snip, "自测·应漂移（片段在别行）"),
        ("bn_core.py", hit_line, "这段代码根本不存在于任何文件 xyzzy", "自测·应未命中"),
        ("bn_no_such_file.py", 3, "print(1)", "自测·文件不存在"),
    ]
    rows, stats = run(samples)
    print(f"\n[1] 分类统计 = {stats}")
    ok &= (stats["命中"] == 1 and stats["行号漂移"] == 1 and stats["未命中"] == 1 and stats["文件不存在"] == 1)
    ok &= all(r[3] is not None for r in rows if r[2] == "行号漂移")
    print(f"[2] 漂移条已报出实际行号: {all(r[3] is not None for r in rows if r[2] == '行号漂移')}")
    print(f"[3] 未命中条目已显式打印（上面可见）: True")
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    if len(sys.argv) > 1:
        text = open(sys.argv[1], encoding="utf-8").read()
        _, stats = run(parse(text))
        print(f"\n汇总：{stats}")
        bad = stats["未命中"] + stats["文件不存在"]
        print(f"不可采信条数（未命中+文件不存在）= {bad}")
        sys.exit(1 if bad else 0)
    sys.exit(0 if selftest() else 1)

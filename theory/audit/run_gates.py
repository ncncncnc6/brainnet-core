#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""theory/audit/run_gates.py —— 可信闸门与回执（gates/1.0.0）

修测量工具本身：此前 shell 循环写成 `timeout N python3 f.py; echo "$(basename $f) exit=$?"`
——**命令替换会先执行并覆盖 `$?`**，退出码被读成 basename 的 0。
本件用 subprocess 直接取 returncode，并把结果写进带时间戳的回执，支持与上次**差分**。

用法：python3 theory/audit/run_gates.py [--strict] [--only 子串]
判据：[1] 退出码取自 subprocess.returncode（不经 shell）
      [2] 回执含每件 exit/耗时/stdout 尾行
      [3] 差分能在删掉回执条目时报出该件（负例）
"""
VERSION = "gates/1.0.0"
import json, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RECEIPT = os.path.join(ROOT, "theory", "audit", "gate_receipt.json")
SKIP = {"run_gates.py", "verify_claims.py"}

def py_files():
    out = []
    # [自纠] 原只扫 ROOT 与 theory/ → theory/audit/ 下的件（探针/诊断）从未进闸门。
    for d in (ROOT, os.path.join(ROOT, "theory"), os.path.join(ROOT, "theory", "audit")):
        if not os.path.isdir(d): continue
        for f in sorted(os.listdir(d)):
            if f.endswith(".py") and f not in SKIP:
                out.append(os.path.join(d, f))
    return out

def run_one(path, timeout=280):
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, "-u", path], cwd=ROOT, timeout=timeout,
                           capture_output=True, text=True)
        rc, out, err = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        rc, out, err = "TIMEOUT", "", ""
    lines = [l for l in out.strip().splitlines() if l.strip()]
    return {"exit": rc, "sec": round(time.time() - t0, 1),
            "tail": lines[-1][:100] if lines else "",
            "err": err.strip().splitlines()[-1][:100] if err.strip() else ""}

def compile_all():
    bad = []
    for p in py_files():
        r = subprocess.run([sys.executable, "-c",
                            "import py_compile;py_compile.compile(%r,doraise=True)" % p],
                           cwd=ROOT, capture_output=True, text=True)
        if r.returncode != 0: bad.append(os.path.relpath(p, ROOT))
    return bad

def diff(prev, cur):
    added_red, fixed = [], []
    for k, v in cur.items():
        if k not in prev: continue
        if prev[k] == 0 and v != 0: added_red.append((k, v))
        elif prev[k] != 0 and v == 0: fixed.append(k)
    return added_red, fixed

def main():
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    strict = "--strict" in sys.argv
    print(VERSION + "\n")
    comp = compile_all()
    print("[闸门① 编译] 失败 %d 件 %s" % (len(comp), comp if comp else ""))
    files = [f for f in py_files() if (only is None or only in f)]
    res, t0 = {}, time.time()
    for p in files:
        rel = os.path.relpath(p, ROOT)
        r = run_one(p)
        res[rel] = r["exit"]
        flag = "绿" if r["exit"] == 0 else ("红" if r["exit"] == 1 else str(r["exit"]))
        print("  %-4s %-34s %6.1fs  %s" % (flag, rel, r["sec"], r["tail"][:52]))
    reds = [k for k, v in res.items() if v != 0]
    print("\n[闸门② 全量] 绿 %d / 共 %d；红: %s" % (len(res) - len(reds), len(res), reds))
    probes = {}
    for pr in ("theory/audit/probe_integrity.py", "theory/audit/probe_traceability.py"):
        rp = run_one(os.path.join(ROOT, pr))
        probes[pr] = rp["exit"]
        print("[闸门③ %s] exit=%s" % (os.path.basename(pr), rp["exit"]))
    prev = {}
    if os.path.exists(RECEIPT):
        try: prev = json.load(open(RECEIPT, encoding="utf-8")).get("results", {})
        except Exception: prev = {}
    added_red, fixed = diff(prev, res)
    print("\n[差分 vs 上次回执] 新增红 %s；转绿 %s" % (added_red or "无", fixed or "无"))
    # [自纠] 原来 `--only` 也整份覆盖回执 → 过滤跑会把基线写空（下一次差分就永远"无变化"）。
    # 现在：过滤跑只**合并**进已有回执；全量跑才重写。
    merged = dict(prev)
    merged.update(res)
    json.dump({"version": VERSION, "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
               "compile_failed": comp, "results": merged if only else res, "probes": probes,
               "partial": bool(only), "prev_added_red": added_red, "prev_fixed": fixed},
              open(RECEIPT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("[回执] 已写 %s（耗时 %.1fs）" % (os.path.relpath(RECEIPT, ROOT), time.time() - t0))
    ok = (not comp) and (not added_red)
    if strict: ok = ok and (not reds)
    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)

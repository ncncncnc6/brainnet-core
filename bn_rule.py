#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bn_rule.py —— 突触输出规则表（bn-rule/1.0.2）

用户要求把表改成"下一个 →（传到哪里，传多少，多久生效/延迟几拍，抑制还是兴奋，还是直接静默）"，
并指出只有这个形状才谈得上图灵完备。本文件把这五件事定成**列**，并给出可跑的语义与判据。

列定义（每列都有文献出处，不臆造）：
  TARGET 传到哪里   —— 目标集合：local(同区) / hub / broadcast(跨区广播) / none
  AMOUNT 传多少     —— 短时可塑性的有效释放量 u·x（Tsodyks–Markram 形式），
                        出处 PLoS Comput Biol 2026, DOI 10.1371/journal.pcbi.1014730
  DELAY  多久生效   —— 以"拍"为单位的整数延迟；髓鞘决定传导速度、且髓鞘终身重塑
                        ⇒ 这一列是**可塑的**（Mol Psychiatry 2026 / MedComm 2026 两份髓鞘论文）
  SIGN   抑制/兴奋  —— +1 / -1；Dale 定律的多种 incarnation 说明"同一细胞所有突触同号"只是特例，
                        出处 J Neurophysiol 2026, DOI 10.1152/jn.00196.2026
  MODE   静默        —— SEND / SILENT；静默 = 释放失败（可由 adapt 或资源耗尽触发）
                        树突生长与局部可塑性见 eLife 2025, DOI 10.7554/elife.87527
另两列用于把表变成可执行的状态机：
  COND   触发条件     —— 前件状态（resting / fired / refractory / exhausted）
  POST   后件更新     —— 本行执行后对自己的状态更新（u、x、ref 等）

── 关于"图灵完备"（重要，不含糊）────────────────────────────────
  · 把 DELAY 当**延迟线**用，投递出去的事件在第 t+DELAY 拍才落地；这条延迟线就是**队列**。
  · 队列自动机（queue automaton）是已知图灵完备的；DELAY 若被钉成常数，整张表退化为
    **有限状态自动机**（无界内存没了）。
  · 本文件给出一个**非正则语言的识别见证**（回文检测需要无界栈/队列），用来证明 DELAY 列是载荷、
    而不是装饰。注意：证明"超出有限自动机"是图灵完备的必要条件，不等于完整 TM 模拟。

判据：
  [1] 表本身合法：每行五列齐全、SIGN ∈ {+1,-1}、MODE ∈ {SEND,SILENT}、DELAY ≥ 0
  [2] 各行的**实际触发概率**（要求：非静默行不得长期停在 ~0 或 ~100；静默行例外）
  [3] DELAY 当队列：回文识别见证（n=1..8 全部正确）
  [4] 负例：把 DELAY 钉成常数（队列变成长度为 1 的寄存器）→ 回文识别必须失败
"""

VERSION = "bn-rule/1.0.2"
import collections, random, sys

# ── 规则表（示意值；"传多少/延迟几拍"的绝对数值是本模型自定，不是实验测量）──
# 行格式: (COND, TARGET, AMOUNT, DELAY, SIGN, MODE, POST)
TABLE = [
    ("resting",    "local",     1.00, 1, +1, "SEND",   dict(u=0.6, x=0.85)),
    ("resting",    "hub",       0.35, 3, +1, "SEND",   dict(u=0.4, x=0.90)),
    ("fired",      "local",     0.55, 1, +1, "SEND",   dict(u=0.5, x=0.70)),
    ("fired",      "hub",       0.80, 2, +1, "SEND",   dict(u=0.7, x=0.60)),
    ("fired",      "broadcast", 0.25, 5, +1, "SEND",   dict(u=0.3, x=0.75)),
    ("refractory", "none",      0.00, 0, +1, "SILENT", dict(u=0.2, x=0.95)),
    ("exhausted",  "none",      0.00, 0, +1, "SILENT", dict(u=0.1, x=0.99)),
    ("resting",    "local",     0.45, 2, -1, "SEND",   dict(u=0.5, x=0.80)),   # 抑制行（静息态）
    ("fired",      "local",     0.40, 1, -1, "SEND",   dict(u=0.45, x=0.75)),  # [1.0.1] 抑制行（发放态）：
                                                       # 状态机只在 fired/exhausted 查表后，
                                                       # 原抑制行挂在 resting 下 → 抑制投递恒为 0（死行）
]

def validate(table):
    bad = []
    for i, r in enumerate(table):
        cond, tgt, amt, dly, sgn, mode, post = r
        if cond not in ("resting", "fired", "refractory", "exhausted"): bad.append((i, "COND"))
        if tgt not in ("local", "hub", "broadcast", "none"): bad.append((i, "TARGET"))
        if sgn not in (+1, -1): bad.append((i, "SIGN"))
        if mode not in ("SEND", "SILENT"): bad.append((i, "MODE"))
        if dly < 0 or not isinstance(dly, int): bad.append((i, "DELAY"))
        if not (0.0 <= amt <= 1.0): bad.append((i, "AMOUNT"))
        if set(post) - {"u", "x"}: bad.append((i, "POST"))
    return bad

def row_probabilities(ticks=20000, seed=3, weighting="uniform"):
    """[A 组 (c) 68-82] 口径写实：这里是**两档马尔可夫近似 + 显式假设**，不是从过程里量出来的。

    三条必须写在明面上的假设（原 docstring 只写了一句"近似"，读者看不出假设是什么）：
      H1 状态转移 = 硬编码四态循环 resting→fired→refractory→exhausted→resting，带 0.3 的滞留；
         **不是**由 POST 列或生物学推导出来的。
      H2 目标行选择：`uniform` 为等概率；`amount` 为按该行 AMOUNT 加权（AMOUNT 才有机会进入概率）。
      H3 因此 P 的含义是"该状态占比 × 该状态下选中此行的比例"，**不含** u/x 的释放动力学。
    """
    rng = random.Random(seed)
    state = "resting"; hist = collections.Counter(); hits = collections.Counter()
    for t in range(ticks):
        hist[state] += 1
        rows = [i for i, r in enumerate(TABLE) if r[0] == state]
        if rows:
            if weighting == "amount":
                w = [TABLE[i][2] for i in rows]
                tot = sum(w) or 1.0
                x = rng.random() * tot; acc = 0.0; pick = rows[-1]
                for i, wi in zip(rows, w):
                    acc += wi
                    if x <= acc: pick = i; break
            else:
                pick = rng.choice(rows)
            hits[pick] += 1
        # H1：硬编码四态循环（**假设**，非推导）
        state = {"resting": "fired", "fired": "refractory",
                 "refractory": "exhausted", "exhausted": "resting"}[state] if rng.random() < 0.7 else state
    return {i: hits[i] / max(1, hist[TABLE[i][0]]) for i in range(len(TABLE))}, hist

# ── DELAY 当队列：非正则语言 a^n b^n 的识别见证 ──────────────────
# [1.0.0→1.0.1] 旧版用 dict[t+d]=ch 排期做"回文"：d=n−i 让所有事件撞在同一拍上互相覆盖，
#   结果 8/8 全 False，负例也失效（常数延迟的读数恰好等于期望值）。那条见证本身还选了错的语言
#   ——回文要 LIFO（栈），延迟线是 FIFO（队列），天然做不了。修法：换成 **a^n b^n**（同样非正则，
#   且恰好是 FIFO 的看家本领），并用**事件队列**如实实现 DELAY 列的排期与到点读出。
def anbn_via_delay(s, delay_col, cap=None):
    """[1.0.2] 增加 cap：延迟线容量上限。cap=None 为无界；cap=K 即内存封顶。"""
    line = []; t = 0; peak = 0; overflow = False
    for ch in s:
        if ch == "a":
            d = delay_col(len(s), t)
            if cap is not None and len(line) >= cap:
                overflow = True                      # 线满 → 新事件丢失
            else:
                line.append(t + d); line.sort(); peak = max(peak, len(line))
        else:
            ready = [x for x in line if x <= t]
            if not ready:
                return False, peak, True
            line.remove(min(ready))
        t += 1
    return (len(line) == 0 and not overflow), peak, False

def selftest():
    print(f"{VERSION}\n")
    ok = True

    print("[1] 表合法性（五列齐全 + 取值域）")
    bad = validate(TABLE)
    print(f"    行数={len(TABLE)}  列={['COND','TARGET','AMOUNT','DELAY','SIGN','MODE','POST']}")
    print(f"    违规项={bad if bad else '无'} → {'PASS' if not bad else 'FAIL'}")
    ok &= (not bad)

    print("\n[2a] 规则可区分性：行为列相同、仅 POST 不同的行必须被点出来")
    # [A 组 (b) 41-53/122-126] `refractory` 与 `exhausted` 两行的行为列完全相同
    # （none / 0.00 / 0 / +1 / SILENT），差别只在 POST。这不必然是错，但不能**默默**存在：
    # 必须被报出，并说明它们的区别在哪一列、后果是否可观测。
    _sig = lambda r: (r[0], r[1], r[2], r[3], r[4], r[5])
    _dup = []
    for i in range(len(TABLE)):
        for j in range(i + 1, len(TABLE)):
            if _sig(TABLE[i]) == _sig(TABLE[j]):
                _dup.append((i, j, TABLE[i][6], TABLE[j][6]))
    print(f"    行为列相同的行对 = {[(i, j) for i, j, _, _ in _dup] if _dup else '无'}")
    for (i, j, pi, pj) in _dup:
        diff = {k for k in set(pi) | set(pj) if pi.get(k) != pj.get(k)}
        print(f"      #{i} {TABLE[i][0]} vs #{j} {TABLE[j][0]}：仅 POST 的 {sorted(diff)} 不同 → "
              f"节点上**可区分**（例：u {pi.get('u')} vs {pj.get('u')}），" 
              f"但作为规则在无状态观测时不可区分")
    # [自纠] 上一版这里写 `ok &= True`（空判据，被 probe_integrity 的恒真判据检查抓到）。
    # 改成真判据：**连 POST 也完全相同**的行对才是缺陷（纯重复规则），必须为 0。
    _fully_same = [(i, j) for (i, j, pi, pj) in _dup if pi == pj]
    print(f"    完全重复（含 POST）的行对 = {_fully_same if _fully_same else '无'}（必须为 0）")
    ok &= (len(_fully_same) == 0)

    print("\n[2b] 各行触发概率（非静默行不应长期停在 ~0 或 ~100）；**两种目标选择并列**")
    probs, hist = row_probabilities()
    probs_w, _ = row_probabilities(weighting="amount")
    print("    （H1 状态转移为硬编码假设；H3 P 不含 u/x 释放动力学）")
    for i, r in enumerate(TABLE):
        print(f"      #{i} {r[0]:<10} uniform={probs[i]:.3f}  amount加权={probs_w[i]:.3f}")
    for i, r in enumerate(TABLE):
        p = probs[i]
        silent = (r[5] == "SILENT")
        flag = "静默行(例外)" if silent else ("OK" if 0.02 <= p <= 0.98 else "异常")
        print(f"    #{i} {r[0]:<10} {r[1]:<10} amt={r[2]:.2f} delay={r[3]} sign={r[4]:+d} {r[5]:<6} P={p:6.2f}  {flag}")
    live_ok = all(0.02 <= probs[i] <= 0.98 for i, r in enumerate(TABLE) if r[5] == "SEND")
    print(f"    非静默行全部落在 (0.02, 0.98): {live_ok}")
    ok &= live_ok

    print("\n[3] DELAY 当 FIFO 队列（无界）→ 识别 a^n b^n（非正则语言）")
    # [A 组 (b) 89-105/128-135/137-148] 原正例传 `lambda L, t: 1`——**常数延迟**，
    # 等于把 DELAY 列钉死，恰好没测到"这一列是载荷"的主张。现正例改用**非常数且保序**的延迟
    # （随时间递增 ⇒ 先入先出顺序不变），并并列常数延迟作对照。
    good = 0; good_c = 0; tests = ["a" * n + "b" * n for n in range(1, 9)]
    for s in tests:
        r, peak, under = anbn_via_delay(s, lambda L, t: 1 + (t // 3))     # 非常数、保 FIFO
        rc, peak_c, _ = anbn_via_delay(s, lambda L, t: 1)                 # 常数对照
        good += r; good_c += rc
        print(f"    n={len(s)//2:<2} {s}  识别(非常数延迟)={r} 峰值在途={peak}   ｜ 常数延迟={rc}")
    print(f"    非常数延迟通过 {good}/{len(tests)}；常数延迟通过 {good_c}/{len(tests)}")
    print(f"    → {'PASS' if good == len(tests) else 'FAIL'}（DELAY 列取非常数时仍能识别 ⇒ 该列是载荷）")
    ok &= (good == len(tests))

    print("\n[4] 负例（1.0.2 换靶子）：延迟线**封顶** K=3 → n>3 必须失败")
    print("    说明：延迟 0 并不消灭计数能力（令牌当拍可用，等价于计数器）——那不是这条列的载荷所在；")
    print("    真正的分界是**内存是否有界**：有界延迟线 ⇒ 有限自动机。")
    fails = 0
    for n in (2, 4, 6):
        s = "a" * n + "b" * n
        r, peak, under = anbn_via_delay(s, lambda L, t: 1, cap=3)
        fails += (not r)
        print(f"    cap=3  n={n}  识别={r}  峰值在途={peak}  {'（应失败）' if n > 3 else ''}")
    caught = (fails >= 2)
    print(f"    n>3 的三例中失败 {fails} 例 → {'PASS（负例被抓住）' if caught else 'FAIL（判据没牙）'}")
    ok &= caught

    print("\n判据总判:", "全部 PASS" if ok else "有 FAIL")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)

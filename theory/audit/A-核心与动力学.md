# A-核心与动力学 审计报告

审计对象：`bn_core.py` `bn_adex.py` `bn_brown.py` `bn_silence.py` `bn_wave.py` `bn_transmit.py` `bn_rule.py` `bn_units.py`
模式：只读审计（本报告是唯一新建文件，未改动任何源码）
依据：每个文件逐行读全文 + 实跑记录（退出码与关键输出见第 0 节）

## 0. 实跑记录

命令：`cd /sdcard/dsphworkspace/brainnet-core && timeout 200 python3 -u <file>`

| 文件 | 退出码 | stdout 关键行 | stderr |
|---|---|---|---|
| bn_core.py | 0 | `RS: ... 稳态 rate=0.0468（偏差 6.5%）` / `FS: ... rate=0.1363（偏差 2.7%）` / `负例(原配方)：地板=0.0588 < 目标=0.14 → rate=0.0673（偏差 52%）th 漂移=9.998` / `自检总判: 全部 PASS`；放电 186 次，th 3.0 → −3.000 | 空 |
| bn_adex.py | 0 | `dt=0.2/0.1/0.05 脉冲数=11/11/11` / `三者一致=True` / `I_extra=58.11 pA → 实测 20.00 Hz` / `自检总判: 全部 PASS` | 空 |
| bn_brown.py | 0 | `实测指数 = 0.515` / `P(max≥10)=0.4540 2·P(B>=10)=0.4765 差=0.0225` / `μ=0.20 L1 距离=1.4218 → PASS（判据有牙）` / `判据总判: 全部 PASS` | 空 |
| bn_silence.py | **1** | `SILENT_MAX=2.0 发放=186 静默=0` / `SILENT_MAX=0.15 发放=177 静默=145 静默概率=45.03%` / `k=8 拍 静默比例=50.30%` / `k=16 拍 静默比例=0.62%` / `S3 落在中间的窗口长度: []` / `判据总判: 有 FAIL` | 空 |
| bn_wave.py | 0 | `Fig4a 厚度中位数=30.0574 pA` / `Fig4h_5 厚度中位数=48.0827 pA` / `每圈周期=48/128/240 拍 → PASS` / `自检总判: A/B 两个命题通过` | 空 |
| bn_transmit.py | 0 | `发出计数 = [16,16,16,16,7,7,7,7,13,13,13,13]` / `交付=142 静默=56 静默占比=0.28` / `树突半径=2.20 mm 覆盖比例=90.7%` / `单向拓扑下从未收到的节点 = [6]` / `判据总判: 全部 PASS` | 空 |
| bn_rule.py | 0 | `非静默行全部落在 (0.02, 0.98): True` / `通过 8/8` / `cap=3 n=4 识别=False, n=6 识别=False` / `判据总判: 全部 PASS` | 空 |
| bn_units.py | 0 | `适应比 1.5524` / `发放率地板 1/(RECOVERY+1) = 0.0588` / `→ 当前参数相容: False 缺口倍数 = 2.605×` / `RECOVERY 16 → ≤6（发放率地板升到 0.1429，可用）` | 空 |

补充实证命令（均在 `/sdcard/dsphworkspace/brainnet-core` 下执行）：

- 跨目录：`cd /tmp && python3 -u /sdcard/dsphworkspace/brainnet-core/bn_silence.py` → `FileNotFoundError: [Errno 2] No such file or directory: '/tmp/bn_core.py'`，退出码 1；`bn_wave.py` 同样报 `/tmp/bn_adex.py`。
- `th` 漂移终值：`BROKEN` 配方跑 20000 拍后 `th=-37.0`，再跑 20000 拍后 `th=-77.0`。
- 高度稀疏发放算例：`bn_adex.run(Fig4h_5, I_extra=3000)` 得到 `w=2763.467`，无任何报错或告警。

## bn_core.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 86-90 | 运行期缺陷 | `sustained` 分配 `x=[0.0]*(horizon+d+2)`（horizon=4000 时为 4004 个），但循环只填到 `x[horizon-1]`（索引 0..3999）；`tail` 取 `x[horizon-100:horizon]` 落在已填区，`head` 取 `x[1:50]` 也落在已填区 | 现代码头部/尾部都在初始化区内，未初始化的 0 未参与；但整个区间右端 `x[4000:]` 被分配后从未使用，属资源与意图不一致的隐患：任何把 `tail` 右移到 `horizon+d+2` 的改动都会立刻读到大片 0 并静默改变判据结果 | `python3 -c "import bn_core as c; print(len([0.0]*(4000+5+2)), c.sustained(0.82,0.2,5))"` | 低 |
| 80-81 | 运行期缺陷 | `if abs(self.rate - p["TARGET_RATE"]) > 1e-6: self.th += TH_STEP*(...)`，`th` 无任何上界/下界/告警。实测负例配方 `RECOVERY=16,TARGET=0.14`：20000 拍后 `th=-37.0`，再 20000 拍后 `th=-77.0`（约 −0.002/拍线性漂移） | 目标不可达时 `th` 无界漂移且不报错，`eff` 随之变成大负数（实测 `eff=-38.821`），门控/阈值语义被静默破坏；调用方看到的是「跑完了、rate=0.0673」这类正常返回值，拿不到「此配置物理不可达」的信号 | `python3 -c "import bn_core as c; p=dict(c.P); p.update(c.BROKEN); n=c.Neuron(p); [n.step(3.5,True) for _ in range(20000)]; print('th',n.th); [n.step(3.5,True) for _ in range(20000)]; print('th',n.th)"` | 高 |
| 109-116 | 判据缺陷 | 首轮豁免是**整圈**豁免：`if first_lap_exempt and lap == 0: continue` 跳过了整圈的门控检查，只在第二圈开始检查。实测 `ring_survives(6,8)=(True,96)`、`ring_survives(6,16)=(True,192)`，手查 `gate_open(8)=True`、`gate_open(56)=True` | 注释与留痕（第 17 行、第 183 行）声称"缺首轮豁免"被修正、并给了"断点位置断言"，但豁免粒度是圈而非首个拍点。对 d 恰为 `PACE_PERIOD` 倍数的参数，第二圈到达点仍全部落在门开拍上，于是"整圈豁免"永远不会被任何算例暴露 | `python3 -c "import bn_core as c; print(c.ring_survives(6,8), c.ring_survives(6,16), c.gate_open(8), c.gate_open(56))"` | 中 |
| 119 | 运行期缺陷 | `I = math.log2(n_classes)` 被计算出来但从未被使用；`容量格` 直接用 `N <= 2 ** I`，而 `2**log2(16)==16` 恒等 | `I` 是死变量；容量格写成 `N <= 2**I` 使判据与 `n_classes` 参数**恒等**，`n_classes` 除重写表达式外不影响任何结果，看似参数化实为硬编码 | `python3 -c "import bn_core as c; print([x[1] for x in c.five_gates(n_classes=16)], [x[1] for x in c.five_gates(n_classes=4)])"` | 低 |
| 122-123 | 判据缺陷 | `W=24`、`1/df=240` 时窗口格 `W >= 1/df` 判定为 `24 >= 240` → 恒 False；`N=16`、`0.15*W=3.60` 时通道格 `N <= 0.15*W` → 恒 False | 两条 `FAIL` 与任何输入无关，是写死的常数比较；自检第 197 行把这两个 `False` 当作"预期 FAIL 齐备"的正面证据（`neg_ok = (d_["通道格"] is False) and (d_["窗口格"] is False)`），于是"负例"变成恒真断言，判据失去区分能力 | `python3 -c "import bn_core as c; print(c.five_gates())"` | 中 |
| 150-165 | 文件头文字与实际行为不一致 | 第 9-21 行留痕声称 §61"纯延迟周期"与 §1.2"门控存活"已拆成两个函数且"五格口径"逐条修过；实际 `ring_survives` 的豁免粒度（见上）与留痕描述不符，第 197 行的负例是常数恒真 | 留痕是审阅者判断"这一版修好了"的唯一依据；当留痕断言与代码实际行为不一致时，代码审计结论会被留痕误导 | `grep -n "首轮豁免" bn_core.py` 对照 `sed -n '109,116p' bn_core.py` | 中 |
| 132-137 | 运行期缺陷 | `_run(neuron, I, ticks, gate=True)` 的 `gate` 恒为布尔真，`Neuron.step(I, gate_open, pace_gate=None)` 的第 2 位形参名为 `gate_open`，但函数体从未读它（第 70 行只读 `pg`），第 59 行 `pg = NET["PACE_GATE"] if pace_gate is None else pace_gate` 读的是模块级 `NET` | 调用方传入的门控状态被静默忽略：`_run` 里"关门"的写法（传 `False`）与"开门"完全等价，任何依赖 `_run(..., gate=False)` 做对照实验的调用都会得到错误结果而不报错 | `python3 -c "import bn_core as c; n=c.Neuron(); print(len(c._run(n,3.5,300,gate=False)), len(c._run(c.Neuron(),3.5,300,gate=True)))"` | 中 |
| 190 | 判据缺陷 | `t_ok = (t == n_ * d + d) if not exp else True`：`exp=True`（存活）分支被强制写成 `True` | `ring_survives` 在存活时是跑满 `laps*n` 步后正常返回，必然 `t = laps*n*d` 而非 `n*d+d`；也就是说该断言在能通过的一侧**结构上不可达**。留痕第 18 行称"加 t = n·d + d 位置断言"堵住假绿，实际只堵住死亡分支 | `python3 -c "import bn_core as c; print(c.ring_survives(6,8))"` 对照 `sed -n '188,192p' bn_core.py` | 中 |

## bn_adex.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 101-112 | 度量定义缺陷 | 反解判据 `abs(r - target_hz) <= tol` 取 `tol=0.5`，但 `rate_hz` 的分辨率是 `1000/(T-skip) = 1000/(1000-200) = 1.25` Hz（`dt=0.2` 时只可能给出 1.25 Hz 的整数倍）。实测 `solve_I(...,20.0)` 返回 `I_extra=58.11 pA → 20.00 Hz`，恰好命中，但只是 `20.0/1.25=16` 这个格点的运气 | 判据容差（0.5 Hz）小于度量本身的最小步长（1.25 Hz），意味着存在大量目标值（如 10.6 Hz）在**任何**电流下都返回 `abs(r-tol)>tol`，函数悄悄跑满 60 次二分后返回中点，调用方无法区分"解到"与"没解到" | `python3 -c "import bn_adex as a; print(a.solve_I(a.row('Fig8-RS'),10.6))"` 对照 `sed -n '101,112p' bn_adex.py` | 中 |
| 43 | 文件头文字与实际行为不一致 | 第 43 行注：`倒数第二列原文表头写作 C (β)，含义本次未在正文中找到解释，原样保留、不解读`；第 57 行 `V_PEAK = 0.0  # [论文] 重置阈值固定为 0 mV`，但 `TABLE1` 第 10 列（`Cbeta`）在 Fig8 三行是 8.4/10.4/10.4，Fig4 六行是 `None` | `Cbeta` 列被"解释为未解释"地保留，而 `run()` 解包时用 `_cb` 接收后丢弃（第 68 行 `_, C, gL, EL, VT, dT, a, tw, b, Vr, _cb, I0 = p`）——即该列既不参与计算也不参与任何校验，属"表里有一列、语义未定、静默忽略"，读者会误以为它已被使用 | `python3 -c "import bn_adex as a; print(a.TABLE1[0][10], a.TABLE1[3][10])"` | 低 |
| 84-98 | 运行期缺陷 | 指数项 `A = gL*dT*math.exp((v-VT)/dT)` 无上界保护；实测 `run(Fig4h_5, I_extra=3000, T=500, dt=0.2)` 返回 `w=2763.467`（正常 Fig8-RS 末 `w≈38.7`），无任何报错/告警/有限性检查 | `w` 被允许增长到远超模型量级（`b=30 pA`、`a=-11 nS` 的适应变量）而不触发任何信号；代码没有 `isfinite` 守卫，若 v 在某步冲到 `exp` 溢出会静默变 `inf`/`nan` 继续参与积分。本次实测到 2763 仍未溢出（阈值截断在 0 mV 拦住了主路径），但"无守卫"这一事实成立 | `python3 -c "import bn_adex as a; print(a.run(a.row('Fig4h_5'),dt=0.2,T=500.0,I_extra=3000.0)[2])"` | 中 |
| 150-153 | 判据缺陷 | 第 150-151 行 `a0 = run(row("Fig8-RS"), I_extra=-98.0, T=500.0)`；第 153 行打印 `I=0 pA → {len(a0)} 个；I=98 pA → {len(a1)} 个` | 同一份代码里三个口径互斥：TABLE1 的 Fig8-RS 基础电流是 `I=98 pA`，`I_extra=-98` 使实际总电流为 0（这才是静默的原因），而打印文字写成"0 pA / 98 pA"。留痕第 24-29 行反复强调"判据与所测对象对应"，此处恰好是标签与对象不对应 | `python3 -c "import bn_adex as a; r=a.row('Fig8-RS'); print('I0=',r[-1]); print(a.run(r,I_extra=-98.0,T=500.0)[0])"` | 中 |

## bn_brown.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 48-54 | 判据缺陷 | `reflection_error(mu=0.0, sigma=1.0, t=200, a=10.0)` 的默认参数把 `a/ (sigma*sqrt(t)) = 10/sqrt(200) ≈ 0.707` 钉死；自检第 83 行用绝对阈值 `err < 0.03` | 绝对误差阈值在固定 `a,sigma,t` 下才成立。替换驱动刻度（这正是"生长"语境需要做的）后同一判据会失去意义：`sigma` 减半（0.5）时 `2·P(B_t≥a)` 极小、经验频率量化误差占主导，`err<0.03` 变成几乎恒真或恒假。该判据没有任何相对化 | `python3 -c "import bn_brown as b; print(b.reflection_error(sigma=0.5)[2], b.reflection_error(sigma=3.0)[2])"` | 中 |
| 56-70 | 度量定义缺陷 | `bins=12`、`t=400`、`trials=1500`：`hist[min(bins-1, int(frac*bins))]` 把 `frac=1.0` 并进最后一箱；`frac` 的取值格点只有 `401` 个（`sum(v>0)/len(p)`，`len(p)=401`），相邻箱界上的点会因严格大于号 `v > 0` 而整体偏一侧 | 12 箱的理论质量在端点箱是 0.186、中间仅 0.053，样本量最小的箱（端点）只期望 ~279 个；L1 是**未归一化**的箱概率绝对差之和，不同 `bins` 会给出不同 L1（分箱越细 L1 越小）。第 91/99/104 行的阈值 0.12、`> d0*2` 都是对 `bins=12` 而言的经验值，没有任何分箱无关性论证 | `python3 -c "import bn_brown as b; print(b.arcsine_distance(0.0,bins=6)[0], b.arcsine_distance(0.0)[0], b.arcsine_distance(0.0,bins=24)[0])"` | 中 |
| 28-46 | 度量定义缺陷 | `walk`/`rms_displacement`/`scale_exponent` 用 `sigma=1.0`、`mu=0`、`trials=400`、`n1=100,n2=400`；`scale_exponent` 只由两个尺度点算出 `log(b/a)/log(n2/n1)`，实测 `0.515`，判据 `abs(exp-0.5)<0.08` | 两点比值只是斜率的两点估计，`0.08` 的容差把它变成宽松的"像 √t"检查；它无法区分 `t^0.5` 与 `t^0.55`（后者在 100→400 上给出 `0.55`，已在容差外才是运气）。该量不度量"往外扩展"，只是布朗标度，与文件问的"细胞往外扩展"之间没有任何连接件 | `python3 -c "import bn_brown as b; print([b.scale_exponent(seed_s) if False else b.scale_exponent() for _ in range(3)])"` | 低 |

## bn_silence.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 74-82 | 判据缺陷 | `s3_window` 的静默定义是"最近 k 拍内**一次发放都没有**"（`sum(hist[t-k:t]) == 0`）。实测 `k=2 → 87.58%`、`k=4 → 75.16%`、`k=8 → 50.30%`、`k=16 → 0.62%`、`k=32 → 0.00%` | 该判据在 I=10 的固定发放率下只能取到两个极端，"落在 5%–50% 中间区间"的窗口长度实测为空集（`S3 落在中间的窗口长度: []`），直接导致总判 FAIL（exit=1）。要落在中间区间，k 需要取到发放间隔量级以上的连续值，而代码只测 2/4/8/16/32 五个离散值，且定义用的是"零发放"这个刚性事件 | `python3 -c "import bn_silence as s; print([(k, round(s.s3_window(k),4)) for k in (2,4,8,16,32)])"` | 中 |
| 46-53、106-119 | 度量定义缺陷 | `pr = s/(f+s)` 被打印为 `静默概率=45.03%`（`SILENT_MAX=0.15`）等，标题是"S1 饱和静默：触发概率 vs SILENT_MAX"；但 `f` 是发放次数、`s` 是 `silence_events`（每拍被截断且 `adapt≥SILENT_MAX` 时累加），两者是**事件计数**，分母不是拍数 | 列名用"概率"、口径是"被打断次数/（发放+打断）次数"，分母既不是时间也不是试验次数，数值不可解释为任何概率。更关键的是：`s` 的计数方式是"每次被截断都累加"，因此 `SILENT_MAX=0.01` 时单次触发的计数可高达 2190，`s/(f+s)` 被这个计数机制主导，不能跨 `SILENT_MAX` 比较 | `python3 -c "import bn_silence as s; print(s.rate_of(0.01), s.rate_of(0.15))"` | 中 |
| 107 | 判据缺陷 | `ok &= bool(mid_s1) and bool(mid_s3)`：要求 S1 与 S3 都落在 5%–50% 才算通过。实测 `mid_s1=[0.15]`、`mid_s3=[]` | 这是把两个独立判据的**与**作为总判，任何一个为空即整体 FAIL；但 S3 的空集是定义本身造成的（见上条），不是"发现异常"。失败信号因此指向"窗口长度选取"，与文件真正要回答的问题（静默概率是否太低）之间的联系没有在代码里体现 | `timeout 200 python3 -u bn_silence.py; echo $?`（实测 1） | 中 |
| 114 | 运行期缺陷 | `allfire = (f0 == 0) or (sec := rate_of(0.002))`：`allfire` 被赋值后从未使用；`sec` 的赋值也是死代码 | `sec` 会真实触发第二次 `rate_of(0.002)`（每次调用都重新 `exec_module` 加载 `bn_core`）的副作用被留在自检路径里，既是性能浪费，也是"判据没写完"的痕迹（`allfire` 显然是某个被删掉的断言） | `grep -n "allfire" bn_silence.py`（仅第 114 行一处） | 低 |
| 32-39 | 度量定义缺陷 | `rate_of(sm, I=10.0, ticks=3000, p_extra=None)` 每调用一次就 `core()` 一次（`spec_from_file_location` 重新加载 `bn_core.py`）；`s1_table` 9 档 + 负例 3 次，共 12 次重载；`s3_window` 另建 `Neuron(dict(c.P))` 且 `ticks=5000` | 同一进程内反复 `exec_module` 不会复用已加载模块，属"每次取样都重造实验台"；本文件所有读数都来自 `ticks=3000/5000` 的单条轨迹（`seed` 只用于 `p_negative`，`rate_of`/`s3_window` 完全确定性、只有一个样本），不存在重复试验意义上的"概率" | `python3 -c "import bn_silence as s, time; t=time.time(); [s.rate_of(x) for x in (2.0,0.15)]; print(round(time.time()-t,2),'s')"` | 低 |

## bn_wave.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 61-71、97 | 度量定义缺陷 | `curve_thickness` 返回 `max(w)-min(w)` 的中位数，打印单位写成 `pA`；实测 `Fig4a=30.0574 pA`、`Fig4b=246.8548 pA`、`Fig4h_5=48.0827 pA`。而 `w` 的方程是 `τ_w dw/dt = a(V−E_L) − w`，`a` 的单位是 nS，故 `w` 的数值单位是 `gL·dT` 量纲（nS·mV），不是 pA | 单位标注与量纲不符；更实质的是这个"厚度"对分箱数 `Nbins=40` 和采样步长 `dt=0.2`/`T=300` 强依赖：改变 `Nbins` 或 `dt` 会改变每个箱内的极差，因此 `res["Fig4a"] < res["Fig4h_5"]` 这个比较没有尺度无关性。论文原话说的是"细而连续的曲线"，代码量的是"以 V 分箱后 w 的极差"，两者只在特定分辨率下才近似同一件事 | `python3 -c "import bn_wave as w, bn_adex as a; t,V,W=w.trace_vw(a.row('Fig4a'),dt=0.2,T=300.0); print(w.curve_thickness(V,W,Nbins=10), w.curve_thickness(V,W,Nbins=40), w.curve_thickness(V,W,Nbins=80))"` | 中 |
| 46-58 | 运行期缺陷 | `trace_vw` 在放电分支里 `break`（第 56 行），跳出内层 `while` 循环 | 与 `bn_adex` 留痕第 25-26 行明确修过的缺陷同类：放电后 `break` 会丢掉该外层步剩下的子步时长，导致脉冲数随 `dt` 变化。两个文件是同一套方程的两个实现，`bn_adex/1.1.1` 修了、`bn_wave` 没修，于是本文件的 (V,w) 轨迹在放电拍附近与 `bn_adex` 不一致 | `python3 -c "import bn_wave as w, bn_adex as a; print(len(w.trace_vw(a.row('Fig8-RS'),dt=0.2,T=300.0)[1]), len(a.run(a.row('Fig8-RS'),dt=0.2,T=300.0)[0]))"` | 中 |
| 73-84 | 文件头文字与实际行为不一致 | `ring_wave` 的每圈周期用硬编码 `T0=8`（第 73 行），门控判定 `(t % T0) != 0`；文件头第 13-14 行承认"环、bn_adex 单细胞是给定工程件"，第 15-18 行把命题 C 标为未验证 | 第 104-117 行的结论把"环上行波 = 循环节"写成已成立的等价，但 `T0=8` 是模块级常量、与 `bn_core` 的 `PACE_PERIOD`/`PACE_GATE` 是两套独立约定（`bn_core` 还有 `PACE_GATE=2.0` 的阈移），此处只按拍点开闭复刻。而且该结论跨到了单细胞相平面（`bn_adex`）与网络环（本文件）两个不同载体上，二者没有共享的动力学方程 | `python3 -c "import bn_wave as w; print(w.ring_wave(8,16))"` 对照 `sed -n '104,117p' bn_wave.py` | 低 |

## bn_transmit.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 98-111 | 度量定义缺陷 | `coverage` 内 `if not any(...): continue` 之后紧接着 `tot += 1; hit += 1`，两行无条件成对出现，恒有 `tot == hit`；返回 `hit/total_pts`，其中 `total_pts = samples*samples = 3600`。实测打印 `树突半径=2.20 mm 覆盖比例=90.7%`；总判 `ok &= (cov > 0.3)` | 返回量在数学上等于"图框内采样点数 / 采样总数"，与树突半径、节点位置、边**完全无关**。实测 `reach=min(L_MAX, G_GROW*120)` 对图框内所有点都成立（`X = sph.hypot(x - p[0], y - p[1]) <= reach` 常量真），所以 `any(...)` 这个唯一的树突判据从未否决过任何一点。把 `reach` 改成 0 也会得到同一数字——覆盖比例 90.7% 是伪指标，"空间侧被树突覆盖的面积占比"这一结论无证据支撑 | `python3 -c "import bn_transmit as t, bn_spatial as sp; e,pos=t.edges_of(sp); print(t.coverage(sp,e,pos))"` | 高 |
| 113-133 | 判据缺陷 | `render` 里的行选择是 `rows = [row for row in r.TABLE if row[0] == ("fired" if (t % 8) >= 4 else "resting")]`，即按 `t%8` 硬造状态；而 `simulate`（第 51-96 行）用的是"由节点自身历史推出"的 `state = "exhausted" if ... else "fired"`（第 82 行） | 留痕第 53-55 行明确把"状态由 `t % 8` 硬造"列为 1.0.0→1.0.1 已修缺陷，并说该错"与 `bn_core` 的 `SILENT_MAX` 死分支同一类"。`simulate` 修了，`render` 没修，于是动图与统计用的是两套状态机；动图永远不会呈现 `exhausted` 行，也永远不会呈现 `refractory` 行 | `grep -n "t % 8" bn_transmit.py`（第 126 行一处，`simulate` 内无） | 中 |
| 69-72 | 运行期缺陷 | 到点投递用 `[q for q in queue if q[0] == t]` 全表扫描，随后 `queue = [q for q in queue if q[0] != t]` 再做一次全表扫描；`queue.append((t+dly, ...))` 的延迟 `dly` 取自 `r.TABLE`（0–5 拍） | 当前 `beats=64`、12 节点下无性能问题，但 `delivered`/`in_cnt` 的计数依赖这两次扫描的**顺序**：先统计、后清理。若同一拍有多个同源事件，`by_area` 的累加顺序即影响 `mat` 字典的迭代顺序（第 199 行再累加一次），而第 201-204 行的 `missing` 判定只关心键存在性，故当前不影响结论。属可读性与规模相关的隐患，本次实测未复现错误 | `python3 -c "import bn_transmit as t, bn_spatial as sp; e,pos=t.edges_of(sp); print(t.simulate(e)[2:4])"` | 低 |
| 26-30 | 文件头文字与实际行为不一致 | 第 27-30 行注释给出"实测扫描（每次发放的静默占比）：N=1→100% … **N=8→28.4%**，N=12→0%"并据此取 `EXHAUST_N=8`。实测 `simulate` 输出 `交付=142 静默=56 静默占比=0.28`，与 28.4% 接近；但注释把 `N=1→100%` 这类值写成扫描结果，而扫描代码（`EXHAUST_N` 的七档对照）在本文件中不存在，无法复现该表 | 文件头用一个在本文件内无法复现的扫描表作为"取 8 的理由"，读者无法核对 N=1/2/4/6/12 各档读数是否确实来自当前 `simulate`（`recent` 窗口是 24 拍，注释里写的是"近 24 拍"，第 28 行与第 82 行一致；但第 55 行又写"近 16 拍内发放 ≥4 次"） | `grep -n "16 拍\|24 拍\|EXHAUST_N\|RECOVERY_BEATS" bn_transmit.py` | 低 |

## bn_rule.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 68-82 | 度量定义缺陷 | `row_probabilities` 的 docstring 自述"状态用两档 Markov 走一遍的近似"，转移是硬编码链条 `resting→fired→refractory→exhausted→resting`（第 80-81 行，`rng.random()<0.7` 才前进），目标选择是 `rng.choice(rows)` 均匀抽样 | 所谓"各行实际触发概率"与任何动力学无关：状态占比由马尔可夫链的遍历性决定，行内选择均匀，于是非静默行的 `P` 只由"该状态下有几行"决定（实测 #0/#1/#7 都是 `0.33`、#2/#3/#4/#8 都是 `0.25`）。第 124 行的判据 `0.02 <= p <= 0.98` 因此对非静默行**恒真**（2 行时 0.5、3 行时 0.33，都不可能越界），无法发现"某行永远选不中" | `python3 -c "import bn_rule as r; p,h=r.row_probabilities(); print(p); print(h)"` | 中 |
| 41-53、122-126 | 判据缺陷 | `TABLE` 里 `refractory` 与 `exhausted` 两行都是 `MODE="SILENT"`；第 121-122 行 `silent = (r[5]=="SILENT")`，`flag = "静默行(例外)"`，第 124 行 `if r[5] == "SEND"` 才纳入判据 | 判据显式豁免了静默行。而 `refractory`/`exhausted` 恰是历史上两行"永远选不中"的死行（文件头第 50-52 行留痕承认旧版抑制行挂错状态导致死行），现在它们仍被"静默行例外"豁免，仍不在检测范围。第 5 行第 6 列的 `MODE` 列语义（SEND/SILENT）被用来做豁免理由，等于用待测对象给自己免责 | `python3 -c "import bn_rule as r; print([(i,x[0],x[5]) for i,x in enumerate(r.TABLE)])"` | 中 |
| 89-105、128-135、137-148 | 判据缺陷 | `anbn_via_delay` 的"无界队列识别 a^n b^n"测试（第 131 行）传入 `lambda L, t: 1`——常数延迟 1。所以"延迟线当队列"的载荷测试实际只用了"容量计数器"这一种实例；负例（第 143 行）用 `cap=3` 让延迟线有界，测出 "n>3 失败" 符合"有界内存 ⇒ 有限自动机" | 正例（常数延迟 1）与负例（同一常数延迟 1 加 cap）都不能区分"延迟列真的被当成可变延迟用"与"延迟列只是被读成常数 1"：把 `delay_col` 的实现换成 `return 1`（忽略传入的 `L,t`）会得到逐位相同的结果。文件头第 23-27 行声称"非正则语言的识别见证用来证明 DELAY 列是载荷"，但该见证只证明"容量无界"是载荷，没有证明"DELAY 的数值变化"是载荷 | `python3 -c "import bn_rule as r; print(r.anbn_via_delay('aabb', lambda L,t: 1)); print(r.anbn_via_delay('aabb', lambda L,t: 99))"` | 中 |
| 84-88 | 文件头文字与实际行为不一致 | 第 84 行注释写"非正则语言 **a^n b^n** 的识别见证"，但第 85-87 行紧接着保留旧版说明"旧版用 dict[t+d]=ch 排期做**回文**…那条见证本身还选了错的语言——回文要 LIFO（栈），延迟线是 FIFO（队列），天然做不了" | 留痕正确指出旧靶子错了，但没有说明新靶子的边界：a^n b^n 用 FIFO 识别需要"在 b 阶段按到达顺序弹出"，而 `anbn_via_delay` 用 `min(ready)` 弹出，对常数延迟等于 FIFO、对可变延迟则不保证。这层前提没有写进留痕 | `python3 -c "import bn_rule as r; print(r.anbn_via_delay('aabbb', lambda L,t: t+1))"` | 低 |

## bn_units.py

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
|---|---|---|---|---|---|
| 87、90-91 | 文件头文字与实际行为不一致 | 第 87 行打印 `当前参数相容: False   缺口倍数 = 2.605×`；第 90 行却写"① 缩不应期：RECOVERY 16 → ≤6（发放率地板升到 **0.1429**，可用）"。同一文件第 105 行 `assertions` 判 `abs(max_rate()-1/17)`，即 `max_rate()` 是 **1/17=0.0588235**，不是 0.1429；0.1429 来自 1/(6+1)，对应改后 `RECOVERY=6` 的新地板 | 同一份报告里"地板"这个词被用于两个不同的量（改前 1/17、改后 1/7），而第 90 行的括号把 0.1429 说成"升到"，读者会把它当成当前地板。第 93 行"①与②同时做会过修（地板 **0.0588** 已 ≥ 目标时…）"又把 0.0588 称作地板。三处数字（0.0588 / 0.1429 / 0.0588）与当前 `RECOVERY` 的对应关系互相矛盾 | `python3 -c "import bn_units as u; print(u.max_rate(), 1/7, u.recovery_compatible(0.14))"` | 中 |
| 21-22、25、78 | 文件头文字与实际行为不一致 | 第 23 行结论写"当前 TARGET_RATE=0.14、RECOVERY=16：16 > 6.142 ⟹ 不相容（**差 2.38 倍**）"；实测 `P['RECOVERY']/ub = 16/(1/0.14-1) = 2.6047`，第 78 行自己打印 `缺口倍数 = 2.605×`。同一文件同一事实给出 2.38 与 2.605 两个数 | 2.38 无来源（16/6.142 的倒数约为 2.605）。这个数被写进文件头作为概览结论，与正文实测值不一致，属可直接核对的算术不一致 | `python3 -c "import bn_units as u; print(16/u.recovery_compatible(0.14), 16/6.142)"` | 低 |
| 42、56-58 | 运行期缺陷 | `bn_units` 内置自己的 `P = dict(DECAY=0.82, U_DECAY=0.88, U_COUPLE=1.0, RECOVERY=16, TARGET_RATE=0.14)`；而同期 `bn_core.PROFILES` 已是 `{"RS": {"RECOVERY":16,"TARGET_RATE":0.05}, "FS": {"RECOVERY":6,"TARGET_RATE":0.14}}`，且 `bn_core.py` 第 18-21 行留痕说明"原配方 RECOVERY=16 + TARGET=0.14 已被证明物理不可达" | `bn_units` 是给"方针"提供决策数字的文件，但它读的是一份与 `bn_core` 已演进的配方**不共享**的本地副本：`bn_core` 当前默认 `P` 的 `RECOVERY=16/TARGET=0.14` 只是遗留默认（真实方针在 `PROFILES`），而 `bn_units` 用 `--assert` 会把这个遗留组合当成"当前参数"去推方针。两个文件对"当前配方"的定义没有单一真源 | `python3 -c "import bn_units as u, bn_core as c; print(u.P); print(c.PROFILES); print(c.P['RECOVERY'], c.P['TARGET_RATE'])"` | 中 |
| 123 | 运行期缺陷 | `sys.exit(0 if ok else 0)`，注释写"本文件是分析件：不相容≠文件错，故仍 exit=0" | 退出码恒定 0，`report()` 的返回值（实测 `False`）被丢弃。任何用退出码串联的自动化（CI、批量脚本）都无法从本文件得知"当前参数不相容"这一结论；同一组的 `bn_core`/`bn_adex`/`bn_rule`/`bn_transmit`/`bn_wave`/`bn_brown` 都用退出码表达判定结果，本文件是唯一例外且无标记 | `timeout 200 python3 -u bn_units.py; echo $?`（实测 0，尽管打印"相容: False"） | 中 |
| 96 | 文件头文字与实际行为不一致 | 第 96 行打印 `判据行的可跑检查：python3 bn_units.py --assert`；第 121-122 行确实实现了 `--assert` → `assertions()`。该分支是存在的，但正文第 1-37 行的"三条与 dt 无关的判据"并未说明这些判据只在 `--assert` 下执行，默认路径只打印 | 这是文档与入口的对应问题而非错误：默认路径返回的是报告、`--assert` 才跑判据。风险在于 `--assert` 里的 `abs(r-1.5524)>1e-4` 等阈值全部是**对 0.82/0.88 反推**的常量（第 103-106 行），一旦 `P['DECAY']`/`P['U_DECAY']` 变化，`assertions` 会 FAIL 但 `report` 不会，两个入口判据不一致 | `python3 bn_units.py --assert; echo $?; python3 bn_units.py; echo $?`（实测 0 与 0） | 低 |

## 本组小结

**缺陷总数：22 条**（其中 1 条记为「未证实」，见下）。

按类别计数：

| 类别 | 条数 | 文件分布 |
|---|---|---|
| (a) 运行期缺陷 | 6 | bn_core 3、bn_adex 1、bn_silence 1、bn_units 2 |
| (b) 判据缺陷 | 8 | bn_core 3、bn_adex 2、bn_brown 1、bn_silence 1、bn_rule 1 |
| (c) 度量定义缺陷 | 6 | bn_adex 1、bn_brown 1、bn_silence 1、bn_wave 2、bn_transmit 1 |
| (d) 文件头文字与实际行为不一致 | 5 | bn_core 1、bn_adex 1、bn_transmit 1、bn_rule 1、bn_units 2 |
| (e) 依赖当前工作目录 | 2 | bn_silence 1、bn_wave 1（亦见于 bn_transmit，同模式） |

注：部分条目同时命中两类，上表按条目主因归类；因此分类计数（6+8+6+5+2=27）大于条目数 22。

**最严重的三条**

1. `bn_transmit.py:98-111`（度量定义，高）——`coverage()` 返回量恒等于"图框内采样点占比"，与树突半径、节点位置、边无关。实测打印 90.7% 并作为"空间侧被树突覆盖"的证据，总判据此 `cov > 0.3` 通过。把唯一判据 `reach` 改成 0 结果不变，属可直接核对的伪指标。
2. `bn_core.py:80-81`（运行期，高）——目标不可达时 `th` 无界漂移且不报错，实测 20000 拍后 `th=-37.0`、40000 拍后 `th=-77.0`，`eff` 随之下沉为 −38.8，神经元退化为"只要门开就发放"。
3. `bn_wave.py:46-58`（运行期，中）——`trace_vw` 在放电分支 `break`，丢掉该外层步剩余时长；这正是 `bn_adex.py` 留痕第 25-26 行已修过的同类缺陷，同一套方程的两个实现出现了不一致的积分口径。

**本组是否有缺陷曾导致错误科学结论**

判断：**有，且发生在「结论已被写进文件头/自检输出」的层面，而非仅停留在局部读数。**

依据：

- `bn_transmit.py` 实测输出 `树突半径=2.20 mm 覆盖比例=90.7%`，总判打印 `[3] 空间覆盖 …` 后 `ok &= (cov > 0.3)` 通过，文件头第 13 行把"给出修复：保证每个节点都有出边与入边（闭环），再看三个数是否归零"列为交付目标。而 `coverage()` 的 90.7% 与树突覆盖无关（`tot` 与 `hit` 恒等），所以该"覆盖已达标"的结论**已被写入交付物**，但支撑它的量是伪指标。
- `bn_silence.py` 的 S3 判据实测落在中间区间为空（`k=16 → 0.62%`、`k=32 → 0.00%`、`S3 落在中间的窗口长度: []`），文件据此给出"S3 没有中间区间"的**否定结论**并让总判 FAIL。这个否定结论本身与数据一致（我复跑确认了这五个数），但推动它成立的判据定义（"k 拍内零发放"）把结论锁死在"要么 ~90% 要么 0%"两端；因此"静默概率太低"这一回答有一半是判据定义造出来的，不是被测系统的性质。
- `bn_rule.py` 的 `[2] 各行触发概率` 实测打印 `P=0.33/0.25` 等，文件头第 31 行称其为"各行的**实际**触发概率"。这些数来自硬编码马尔可夫链 + 均匀选择，与被测规则表的真实动力学无关；但豁免静默行后判据 `0.02<=p<=0.98` 实测恒真（输出了 `非静默行全部落在 (0.02, 0.98): True`），于是"每行都活得下来、没有死行"的结论被写进自检输出。真实历史（留痕第 50-52 行）恰恰说明死行是该表的主要故障模式。
- 反面证据（说明不是全面性问题）：`bn_core` 的 `BROKEN` 配方被显式判为 FAIL（`地板=0.0588 < 目标=0.14 → … 必须 FAIL: PASS`），`bn_adex` 的 dt 收敛在 0.2/0.1/0.05 三档下实测 11/11/11 一致，`bn_brown` 的负例 μ=0.20 实测 L1=1.4218 确实被抓住。这三处说明"判据没牙"不是本组全局特征，但上述三条中的两条（bn_transmit、bn_rule）确实产出过不成立的正面结论。

**未证实项**

- `bn_wave.py:49` 原拟记为"w 更新用了步首 v、与 bn_adex 不一致"，逐行比对后**不成立**：`bn_adex.py:86` 的 `w_new` 同样使用旧 `v`，两处一致。故该条已从表格中撤下，不作为缺陷。
- `bn_transmit.py` 用 `hash((i, t))` 选目标行，我最初怀疑受 Python 字符串哈希随机化影响；实测两次独立进程给出完全相同的 `静默=56 交付=142` 与逐节点发出计数，且元组内是 int（int 哈希不随机化），**未证实**存在跨进程不可复现，故不计为缺陷。

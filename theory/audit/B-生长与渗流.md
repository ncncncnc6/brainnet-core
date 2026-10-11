# B 组（生长与渗流）代码审计

- 审计对象：`bn_growth.py`、`bn_grow_v2.py`、`bn_grow_time.py`、`bn_grow_rule.py`、`bn_spatial.py`、`bn_perc.py`、`bn_part.py`、`bn_repro.py`
- 工作目录：`/sdcard/dsphworkspace/brainnet-core`
- 审计方式：只读；每个文件通读全文 + 原样实跑 + 跨目录实跑；未改任何被审源码。
- 统一实跑命令：`cd /sdcard/dsphworkspace/brainnet-core && timeout 200 python3 -u <file>`
- 统一跨目录实跑命令：`cd /tmp && timeout 60 python3 -u /sdcard/dsphworkspace/brainnet-core/<file>`
- 类别代码：(a) 运行期缺陷 / (b) 判据缺陷 / (c) 度量定义缺陷 / (d) 文件头与代码不一致 / (e) 依赖当前工作目录。
- 严重度：高 = 会产出错误科学结论或直接改变结论方向；中 = 判据空转、度量口径错、跨目录不可运行；低 = 死代码、措辞/版本不一致、分辨率小瑕疵。

---

## bn_growth.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`存活边数=87  衰减删除=0`、`回路条数=6348`、`对齐比例=0.060`、`零分布：均值=0.084 标准差=0.009`、`生长图 0.060 落在零分布第 0.7 百分位；显著高于随机=False`、`单调不降=True`、`判据总判: 全部 PASS`。
跨目录实跑：退出码 `0`（本文件不导入同级模块）。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 139、146、159 | (b) 判据 | 正向判据 [3] 只算 `better`（L139）、[4] 只算 `mono`（L146），二者都没有并入 `ok`；全文件唯一的 `ok` 更新是 L159 负例 `ok &= (not still)`。实跑 [3] 输出 `显著高于随机=False`，总判仍打印 `全部 PASS`。 | 退出码与 `判据总判` 只反映负例。正向假设（生长图对齐比例显著高于随机）实测为假，仍报全部通过；负例也因同样“不显著”而通过，正负例无法区分，判据整体空转。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_growth.py; echo EXIT=$?`（输出含 `显著高于随机=False` 与 `全部 PASS`、EXIT=0） | 高 |
| 71、84、90 | (b) 判据 | 生长参数 `window=1`（L71），观测滞后 `lag = tv - tu` 只允许 `0 <= lag <= 1`（L84），延迟取 `d % T0`（L90）。于是延迟取值恒为 {0,1}；T0=8，长度 3–5 的回路只有在“全部延迟为 0”时才有 Σd ≡ 0。实跑生长图对齐 0.060，低于同规模随机 0.084。 | 生长机制自身被窗口参数限死为 0/1 延迟，结构上几乎不可能产生 mod 8 对齐；再把“显著高于随机”缺失当作通过，等于把机制失败读成成功。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_growth.py`（对齐比例 0.060 < 零分布均值 0.084） | 高 |
| 113-118 | (b) 判据 | `correct_delays` 每轮只试 ±1，仅当 `best[1] > cur` 才接受，否则恢复 `base`（L118）。 | 逐轮对齐比例“单调不降”是算法接受规则的构造保证，不是被检验的性质；[4] 的 `单调不降=True` 是必真命题，永远无法证伪“逐渐修正是否有效”。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_growth.py`（无论输入如何 `单调不降=True`） | 中 |
| 91-93 | (b) 判据 | `idle=1500`，`T=4000`，节点周期为 2/3/4（L75）且发放持续到 t≈T-1；`lastt` 总是接近 3999，判据 `T - 1 - t > idle` 恒不成立。实跑 `衰减删除=0`。 | 衰减删除机制在任何默认参数下都不可触发，判据 [1] 声称检验“未用衰减”，实际从未执行到删除分支。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_growth.py`（`衰减删除=0`） | 中 |
| 137-141、149-158 | (b) 判据 | “显著”只做单侧：`better = f > mu + 2*sd`（L139）、`still = fb > mu2 + 2*sd2`（L156）。实测 f=0.060 落在零分布第 0.7 百分位（显著更低），仍被当作“不显著更高”。 | 只惩罚“显著更高”，不惩罚方向反转（显著更低）；负例 `fb=0.083`、`mu2=0.083` 恰好相等，通过的原因与正向失败的原因相同，判据没有鉴别力。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_growth.py` | 中 |
| 58-59、95-102、154-158 | (b) 判据 | `aligned_fraction` 在无回路时返回 0.0；`random_control` 也全返回 0.0；则 `sd=0`、`better=False`、`still=False`，`ok` 保持 True。 | 即使“回路数=0”（结构完全没长出来），本文件仍会打印 `全部 PASS`；“没有可测量的回路”这种最坏情况反而不触发任何失败。 | 逻辑复现：令 `cycles_of` 返回空集后跑 `selftest()`；退出码仍为 0 | 低 |

---

## bn_grow_v2.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`λ = 2.000 mm⁻¹`、`v1：连接=132 … 中段占比=0.00`、`v2：连接=34 … 中段占比=0.56`、`偏度=+0.79`、`v1 双峰=True / v2 双峰=False`、`q=0.10：长边占比=0.35 / q=0.00：长边占比=0.12`、`判据总判: 全部 PASS`。
跨目录实跑：退出码 `1`，stderr 末行 `FileNotFoundError: [Errno 2] No such file or directory: '/tmp/bn_spatial.py'`。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 32-34 | (e) 工作目录 | `load(name)` 用 `spec_from_file_location(name, name + ".py")` 裸文件名解析同级模块；`selftest` 在 L105 调 `load("bn_spatial")`。 | 换目录运行即在 `/tmp` 找 `bn_spatial.py`，直接 FileNotFoundError；本组有 4 个文件（另有 bn_grow_rule、bn_perc、bn_repro）同病。 | `cd /tmp && timeout 60 python3 -u /sdcard/dsphworkspace/brainnet-core/bn_grow_v2.py; echo EXIT=$?`（EXIT=1） | 中 |
| 17-19 vs 36-40 | (d) 文件头不符 | 头部 docstring L18-19 声明 λ 由“最近邻间距 0.6 mm 处 P=0.60 ⇒ λ = −ln(0.60)/0.6 ≈ 0.851”定；实际默认 `LAMBDA_LIT = 2.0`（L36）、`LAMBDA = LAMBDA_LIT`（L40）。在 d=0.6 mm 处实际 P=exp(−1.2)=0.301，而非 0.60。 | 文件头给出并解释的参数规则与实际执行值不是同一个；读者按头部复算会得到 0.851，与运行输出的 λ=2.000 对不上。L18-19 未标注已被 L36-40 取代。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -c "import math; print(math.exp(-2.0*0.6))"`（0.3012）对照 stdout `λ = 2.000` | 中 |
| 78-80 vs 82-91 | (c) 度量定义 | 同一份报告里“中段占比”有两套口径：`middle_fraction` 用固定区间 [0.6,1.4]（L78-80）；`bimodal_direct` 用按数据最大值分 10 档后的 `h[4]+h[5]`（L89）。实跑 v2 两者分别给 0.56 与 0.24。 | 同名指标两种定义，且 bin 边界随样本最大值漂移；判据 [1]（>0.10）与判据 [3]（mid<=0.15）实际上在量不同的东西，数字不可互相印证。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_v2.py`（同屏出现 `中段占比=0.56` 与 `中间占比=0.24`） | 中 |
| 29 vs 36、138 | (d) 文件头不符 | `VERSION = "bn-grow-v2/1.0.0"`（L29），但正文注释自称 `[1.0.1]`（L36）与 `[1.0.2]`（L138）的改动留痕。 | 版本号与留痕不一致，无法从 VERSION 判断实际跑的修订级别。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -c "import bn_grow_v2 as m"` 不可用（需同目录直接运行）；直接读 L29/L36/L138 | 低 |
| 93-97 | (a) 运行期 | `dist_model_selection()` 定义后全文无调用点（`selftest` 的 [3] 改成了 `bimodal_direct`）。 | 死代码；且它引用的 `bn_grow_rule.fit_best` 与文件头“复用三模型似然选择”的描述已经脱节。 | `cd /sdcard/dsphworkspace/brainnet-core && grep -n dist_model_selection bn_grow_v2.py`（仅定义行 93） | 低 |
| 20-21 vs 67 | (d) 文件头不符 | 头部称远带“以 q 概率**无条件**连接”，代码是 `if not keep and dd > FAR_BAND and rng.random() < q`（L67），即以“指数抽样失败”为前提。 | 措辞“无条件”与实际条件抽样不符；边际概率虽近似 q，但口径描述错误。 | 直接读 L21 与 L67 | 低 |

专项核查（λ 出处）：头部 L36-38 声称默认 2.0 有出处——PMC12380633 图 5 真值 1.2/1.4/1.6/2.8/3.0/3.2、插图另给 1/2/10 mm⁻¹。已在 `papers/PMC12380633.xml` 命中：图 5 说明含 `values 1.2, 1.4, 1.6, 2.8, 3.0, and 3.2` 与 `(R = 0.93)`；图 1C alt-text 含 `λ_decay = 1 mm−1, 2 mm−1, 10 mm−1`。**结论：λ=2.0 与注释声称的文献出处一致，不是缺陷**；但“取区间内的 2.0 作默认”只覆盖 L18-19 设计值的对照说明未同步，见上表第 2 行。

---

## bn_grow_time.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`λ 设计值=0.851 … λ 尺度反推值=5.117`、`平均度随 λ 单调不增: True`、`长边占比首超 0.5 的反转临界：λ ≈ 5.116855762208991`、`tick=60 回路=2`、`首个回路 tick=60`、`判据总判: 全部 PASS`。
跨目录实跑：退出码 `0`（本文件定义的 `load` 从未调用）。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 8-14 | (d) 文件头不符 | 头部称：本轮把 `papers/` 全库正则扫过 `λ=` / `lambda=` / `decay length`，“**没有任何一条给出显式数值 λ**”，因此“不写具体文献值”。但同库 `papers/PMC12380633.xml` 的图 1C alt-text 明写 `λ_decay = 1 mm−1, 2 mm−1, 10 mm−1`，正文用 `λ_decay = 5 mm⁻¹`，图 5 真值 `1.2/1.4/1.6/2.8/3.0/3.2`。 | 本文件引用的同一篇论文就给出了显式数值 λ；“全库无任何显式数值”与论文库事实相反，直接导致本文件选用 5.117 而 bn_grow_v2 选用 2.0，两文件对同一来源给出互相冲突的默认值。 | `cd /sdcard/dsphworkspace/brainnet-core && grep -o "λ_decay = 1 mm−1, 2 mm−1, 10 mm−1" papers/PMC12380633.xml`；`grep -o "λ_decay = 5 mm" papers/PMC12380633.xml` | 高 |
| 140-144 | (b) 判据 | 判据 [3] 的声明是“中段占比由低到高”，代码 `for tick in (10,30,60,200): … mid, sk = shape(...); print(...)` 只打印，没有 `ok &=`。 | 判据 [3] 无论中段占比升、降还是 NaN，都不影响退出码与“全部 PASS”。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_time.py`（[3] 段无 PASS/FAIL 字样） | 中 |
| 136-138 | (b) 判据 | 判据 [2] 的文字是“边数与**巨连通分量**随 tick 单调不减”，代码 `nondec = all(hist[i][1] <= hist[i+1][1] ...)` 只比较 `hist[i][1]`（边数），巨分量 `hist[i][2]` 未被断言。 | 巨分量的单调性没有检查；实测巨分量 0.00→0.83 后持平，若出现回退也不会被抓住。判据文字与所测对象不对应。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_time.py`（[2] 仅 `边数随 tick 单调不减`） | 中 |
| 32-33 | (c) 度量定义 | 默认 λ 取 `LAM_DERIVED ≈ 5.117`；同一论文形式下 bn_grow_v2 默认取 2.0。实跑本文件 λ=5.117 时平均度 1.17、tick=200 仅 14 条边；λ=2.0 时平均度 2.42。 | 同一文献同一公式，两个文件的默认值相差 2.5 倍；本文件注释称“如实标注”，但对“选哪个值”没有与 bn_grow_v2 对齐，无法复现同一组结论。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_time.py` 对照 `bn_grow_v2.py` | 中 |
| 115-124 | (c) 度量定义 | 扫描网格为 `(0.0, 0.851, 2.0, 5.117, 8.0, 12.0)`（L115），“反转临界”取第一个 `lf>0.5` 的网格点并打印全精度 `λ ≈ 5.116855762208991`（L122-124）。 | 5.116855762208991 只是网格点 5.117 的浮点回显，不是临界值的估计；真实穿越只知落在 (2.0, 5.117) 区间内。假精度会误导读者把网格值当测量值。 | 读 L115、L122；stdout `λ ≈ 5.116855762208991` | 低 |
| 28-30、93-101 | (a) 运行期 | 定义了 `load(name)`（L28-30）但全文无调用点；`shape()` 用 `hi = max(xs)` 归一化分箱（L95），tick=30 时只有 2 条边，分箱边界由该样本最大值决定。 | `load` 是死代码；`shape` 的“中段占比”在小样本下由最大值归一，物理口径随样本漂移（tick=30 只有 2 条边，mid 无统计意义）。 | `cd /sdcard/dsphworkspace/brainnet-core && grep -n "load(" bn_grow_time.py`（仅定义） | 低 |

---

## bn_grow_rule.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`真值=M1 指数衰减 判为=M3 重尾 ✗`、`判对率=0.67`、`实测连接距离做模型选择 → 最优模型：M2 双峰`、`负例 … PASS（判据有牙）`、`判据总判: 全部 PASS`。
跨目录实跑：退出码 `1`，stderr 末行 `FileNotFoundError: [Errno 2] No such file or directory: '/tmp/bn_spatial.py'`。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 80-90 | (c) 度量定义 | `fit_best` 直接比三模型的对数似然取最大（L89），但 M1 只有 1 个自由参数，M2、M3 各有 2 个（L83-84），没有 AIC/BIC 或任何复杂度惩罚。实跑 M1 合成数据（真值 lam=2.0，且 2.0 在搜索网格 L81 内）被判为 M3 重尾，判对率 0.67。 | 参数越多的模型天然更易抬高似然；不惩罚复杂度就无法判断“M2 双峰胜出”是形态证据还是模型容量。判据与“三模型可区分”这一被测命题不对应。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_rule.py`（`真值=M1 指数衰减 判为=M3 重尾 ✗`、判对率 0.67） | 中 |
| 129-135、137-144 | (b) 判据 | 判据 [2]“三个部件各自可跑”与 [3]“对本模型实测距离做模型选择”都只 `print`，没有 `ok &=`。 | 这两条判据永远不会失败；`ok` 只由 L124（判对率≥0.6）与 L153（负例）决定。文件声称 [2][3] 是判据，实际只是报告。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_rule.py`（[2]/[3] 段无 PASS/FAIL） | 中 |
| 122-124 | (b) 判据 | 三分类混淆矩阵的对角命中率阈值硬编码 `rate >= 0.6`；随机瞎猜的期望是 1/3≈0.33。 | 0.6 无来源；实测 0.67 只比阈值高 0.07，而 M1 已被误判，说明该阈值不足以支撑“拟合器可用”的结论。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_rule.py`（判对率 0.67） | 低 |
| 146-154 | (b) 判据 | 负例通过 `MODELS["M2 双峰"] = lambda ...` 临时替换全局模型（L149），再比较两次 `fit_best` 的“模型名不同”（L152），最后在 L151 手工恢复。 | 这不是对“gap 退化”做同一判据的检验，而是“换掉模型定义后名字变了”；且替换/恢复不放在 try/finally，异常时会污染后续调用。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_grow_rule.py`（负例判 PASS 仅凭 `M2 双峰 ≠ M3 重尾`） | 低 |
| 45 | (c) 度量定义 | `m1_exp(d, lam, p0=1.0)` 把幅值 p0 固定为 1.0，`fit_best` 也不搜 p0。 | M1 的值域被钉死在 [exp(−λ·2.2), exp(−0.01λ)]，与 M2/M3 的幅值可调不对称；这进一步帮助了 M2/M3 在似然比较中胜出，属模型对比的隐藏不公平。 | 读 L45 与 L81-84（网格只有 lam） | 低 |

专项核查（配分函数 Z）：**Z 已加**。L63-71 定义 `normalizer(f, params, steps=600)` 对 f 在 [D_LO, D_HI]=[0.01,2.20] 上做中点法积分；L73-74 `ll_norm = loglik − n·log Z`；L86 在 `fit_best` 中调用 `ll_norm`。运行结果不再是“什么都收的 M2 白捡最高似然”，负例也能改口。此项不是缺陷。

---

## bn_spatial.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`t=100 边数=132 最长边=2.158 mm L_MAX=2.2`、`首个回路出现在 t=3；该回路最长边=0.079 mm`、`disc 边=132 回路=22418 对齐比例=0.4979`、`annulus 边=100 回路=1890 对齐比例=0.2254`、`L_MAX=1.8 mm … 跨区 0 → PASS`。
跨目录实跑：退出码 `0`（不导入同级模块）。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 78、134-135 | (b) 判据 | `build` 在 L78 先 `if d > ru or d > L_MAX: continue`，即所有入边都 ≤L_MAX；L134-135 的判据 [1] 再断言 `dmax <= L_MAX + 1e-9`。 | 判据 [1] 由构造保证，永远不可能失败；它没有检验“长度上限是否真的作为机制在起作用”，只是一句废话。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_spatial.py`（[1] 恒 True） | 中 |
| 40-43、78-79 | (b) 判据 | LAMBDA=2.5、P_MIN=0.05 ⇒ 保留条件为 `exp(-d/2.5) >= 0.05`，即 d ≤ 7.49 mm；而本布局最大可能距离是 2.158 mm，`exp(-2.158/2.5)=0.422`。 | 距离衰减判据在全部可达边上都恒为真，从不拒绝任何边；文件头 L20 声称的“距离衰减”在本布局下完全没有作用，实际唯一约束是 L_MAX。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -c "import math; print(math.exp(-2.158/2.5))"`（0.422 > 0.05） | 中 |
| 149-153 | (b) 判据 | 判据 [3]“环带 vs 全盘对照”只有 `print`，没有 `ok &=`。 | 两种模式的边数/回路数/对齐比例怎么变都不影响退出码；“对照”不是判据。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_spatial.py`（[3] 段无 PASS/FAIL） | 中 |
| 97-99、101-107 | (c) 度量定义 | `beat()` 直接 `return 2`，docstring 称“拍 = 从过程量：最快起搏区周期（沿用上一轮口径）”，但函数体没有任何计算；`aligned_frac` 只用该常量做 `s % B == 0`。文件内没有零分布/随机对照。 | “拍”是硬编码常量，却以“从过程量得出”描述；对齐比例因此只是延迟和的奇偶性。disc 实测 0.4979 与随机期望 0.5 几乎重合，无对照就无法判断这是机制对齐还是偶然。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_spatial.py`（disc 对齐比例=0.4979034704255509） | 中 |
| 145 | (b) 判据 | docstring（L27-28）与输出文字要求“上一刻 < 它”（严格大于），代码实现为 `longest > r_prev - 1e-9`（L145）。 | 判据实际放松了 1e-9：当 `longest` 恰等于或微小于上一刻半径时也会判 True，与文字口径不一致。本次实跑 longest=0.079 > r_prev=0.070，尚未触发该放松。 | 读 L145；`cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_spatial.py`（t=3、longest=0.079、r_prev=0.070） | 低 |
| 30 vs 155-167 | (d) 文件头不符 | 头部判据 [4] 写“把 L_MAX 压到小于闭合所需弦长，**回路必须一个都不出现**（exit=1 判 FAIL）”；实际负例（L155-167）改靶为“只需跨区回路=0，区内回路仍可存在”，实跑总回路 42、区内 42。 | 文件头仍保留旧判据口径（全部回路=0），与已实现的新口径（仅跨区=0）不一致；读者按头部会误以为该负例检验全部回路。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_spatial.py`（`总回路=42（区内 42 / 跨区 0）`） | 低 |

补充核验：文件头 L41 称跨区最小弦长 2.041 mm，实测跨区距离 min=2.0406 mm、max=2.1580 mm，L_MAX=1.8 时按 build 规则可建跨区边=0，故该负例的目标设置有效，非缺陷。

---

## bn_perc.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`覆盖比例=62.1%  平均多重度=5.85  最大多重度=12`、渗流曲线 21 档中前 20 档 `巨分量占比=0.33`、`L_MAX=2.10 mm … 巨分量占比=1.00 ← 首次 ≥0.5`、`实测转变点 L_MAX_c = 2.1 mm；ER 判据 ⟨k⟩=1 对应的边密度 = 0.0909`、`[3] 凋亡：**本件不再测量**`、`判据总判: 全部 PASS`。
跨目录实跑：退出码 `1`，stderr 末行 `FileNotFoundError: [Errno 2] No such file or directory: '/tmp/bn_spatial.py'`。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 3-10、17-21 vs 126-130 | (d) 文件头不符 | 头部 L3 标题写“占位 / 渗流阈值 / **凋亡阈值** 三问的量测”，L9-10 写“量位点渗流…给出 f_c（凋亡阈值）”，L20 判据 [3] 写“随机杀节点 0%→90%…给出实测 f_c”；实跑 L126 打印“**本件不再测量**”，L130 仅 `ok &= True`。 | 文件头承诺的量测已在代码中退役但头部未同步；L130 的 `ok &= True` 让判据 [3] 对总判零贡献，条目变成纯占位。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_perc.py`（stdout 同时出现“给出实测 f_c”的头部承诺与“本件不再测量”） | 中 |
| 78-88、91-104 | (a) 运行期 | `apoptosis(edges, n, frac, seed)`（L78-88）与 `sparse_edges(...)`（L91-104）在全文件中没有任何调用点。 | 死代码。凋亡量测“退役”后函数体仍留在文件里；`sparse_edges` 更与“按论文形式建立稀疏连接”的接口重复，容易被后续误调用却无人维护。 | `cd /sdcard/dsphworkspace/brainnet-core && grep -n "apoptosis\|sparse_edges" bn_perc.py`（仅定义行 78、91，无调用） | 低 |
| 127、129 | (d) 文件头不符 | L127/L129 称“N=400 的 P4′ 测量：**f_c ≈ 0.85**”且“唯一出处是 theory/p45_experiments.py（N=400 ⇒ f_c≈0.85）”。实跑 `theory/p45_experiments.py` 输出 `f_c = 0.8`。 | 本文件把唯一出处的数字写错（0.85 vs 0.80）；p45 的网格点恰有 0.8 与 0.85，实跑第一个跌破 0.5 的是 0.8，说明 0.85 是未核对的手抄值。 | `cd /sdcard/dsphworkspace/brainnet-core && timeout 200 python3 -u theory/p45_experiments.py \| grep "f_c ="`（输出 `f_c = 0.8`） | 中 |
| 122-124 | (b) 判据 | 算 `er_pred = 1/(n-1)`（n=12 → 0.0909）并打印“ER 判据 ⟨k⟩=1 对应的边密度”，但代码从未计算实测边密度，也从未与 `er_pred` 做任何比较；判据 [2] 的 `ok` 只断言 `trans is not None`（L124）。 | 判据 [2] 声称“与 ER 判据 ⟨k⟩=1 的预测对照”，实际没有任何对照发生；ER 预测值只是被打印出来，与几何渗流结果（L_MAX_c=2.1 mm）不可比。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_perc.py`（`er_pred` 出现一次后再无比较） | 中 |
| 63-75、115-124 | (c) 度量定义 | `percolation_by_L` 用 `steps=22` 在 [0.10,2.20] 上扫描（步长 0.1）；N=12 时巨分量只能取 0、1/12、…、1。实跑从 0.33 一步跳到 1.00，转变点报 2.1 mm。 | 空间分辨率 0.1 mm + 巨分量量化 1/12，无法定位真正的 L_MAX_c；“转变点 2.1 mm”只说明 2.0 与 2.1 之间发生了跳变，精度不足以支撑论文级阈值。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_perc.py`（2.00→2.10 之间 0.33→1.00） | 中 |
| 34-42 | (c) 度量定义 | `overlap` 在固定窗口 `xs = -3.5 … 3.5`（L35）上取 80×80 采样，返回 `covered/total`；节点实际分布在半径约 1.22 mm 内。 | “覆盖比例 62.1%”把大片与结构无关的空白窗口算进分母，窗口一换比例就变；该指标不是结构的面积归一化量，作为“空间占位”度量口径不稳定。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_perc.py`（`覆盖比例=62.1%`；窗口写在 L35） | 低 |
| 111-113 | (b) 判据 | L113 打印“平均多重度 > 3 即「严重重叠」: True”，但该比较结果没有并入 `ok`；阈值 3 无来源。 | “严重重叠”的判定只是文字，不构成可失败的判据；阈值 3 与 N、半径无关。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_perc.py`（[1] 段无 PASS/FAIL） | 低 |

专项核查（凋亡判据是否退役、死代码是否仍在执行）：**已在 `selftest` 退役**——L126-130 只打印说明，不再调用凋亡测量；`apoptosis()` / `sparse_edges()` 仍在文件里但**无任何调用点，运行时不会执行**（grep 确认仅定义行）。残留问题是：头部仍承诺该判据、`ok &= True` 占位、死函数未删。

---

## bn_part.py

实跑：退出码 `0`；stderr 为空。关键 stdout：`边数=44/132（区内 36、跨区 8）`、`延迟多重集={1: 36, 3: 8} 不同取值=2`、`回路=42（触顶=False） 区内=42 跨区=0`、`对齐比例：区内=0.0 跨区=None`、`intra 实测=0.000 零分布=0.055±0.031 第 0.0 百分位`、`判据总判: 全部 PASS`。
跨目录实跑：退出码 `0`（不导入同级模块）。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 13、18 vs 74、79 | (d) 文件头不符 | 头部 L13 称延迟方案为“区内 1 拍、hub→叶 3 拍、**hub→hub 2 拍**”，L18 称“mod 8 至少跨 **3 个取值**”。但 L74 `if not same and area[u] not in hubs: continue` 只允许从 hub 区出发的跨区边，而 L79 的延迟 2 对应 `area[v] in hubs`（v 也在 hub 区），与 `not same` 互斥——延迟 2 结构上不可达。实跑延迟多重集 `{1: 36, 3: 8}`，只有 2 个取值。 | 头部承诺的三值延迟方案实际只有两值；“至少跨 3 个取值”与运行输出直接矛盾，[B2] 的修复说明被夸大。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_part.py`（`延迟多重集={1: 36, 3: 8}`） | 中 |
| 20、128-143 | (b) 判据 | 头部 L20 称“区内/跨区回路的对齐比例……**各自**与零分布比较”。实跑 `回路=42（区内 42 跨区 0）`，跨区回路恒为 0；`frac_by_class` 给 `f['cross']=None`，L137-143 直接 `continue` 跳过。 | 跨区回路的对齐判据所测对象在本模型输出里不存在（结构性无跨区回路），判据空转；承诺的两类比较只兑现了区内一类，却仍打印“全部 PASS”。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_part.py`（`跨区=0`、`跨区=None`、零分布段只有 intra 一行） | 高 |
| 145-151 | (b) 判据 | 负例声称“把竞争改回跨类混排（重放 B1）”，实际调用 `grow_partitioned(topk_cross=0)`，直接把跨区配额置 0；然后检查跨区边=0（L150）。 | 把要测的量本身置为 0，再用它必须为 0 当判据，是同一断言的重言式；它没有重放 B1 的“跨类混排”实现，无法证明 [B1] 已修复。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_part.py`（`跨区配额=0 时跨区边=0 … PASS`） | 中 |
| 136-143 | (b) 判据 | 零分布段唯一断言是 `ok &= (sd > 0)`（L143）；没有方向或显著性检验。实跑 intra 实测 0.000，在零分布 0.055±0.031 的第 0.0 百分位——即显著**低于**随机，仍判 PASS。 | 判据只查“零分布有没有方差”，不查实测是否落在预期方向；显著低于随机这一反向结果不会被发现。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_part.py`（`实测=0.000 … 第 0.0 百分位` 且全部 PASS） | 中 |
| 116-126 | (b) 判据 | [1]“三个 bug 的修后体检”打印 `跨区非零: True`（L125）、`非退化: True`（L126），但这两项没有并入 `ok`（L136 的 `ok = varied` 只看对齐比例是否非平凡）。 | [B1]、[B2] 的修复结论只以打印呈现，不构成可失败判据；跨区边数归零或延迟退化成单值时，这两行仍会照常打印而退出码不变。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_part.py`（[1] 段无 PASS/FAIL 且不影响 EXIT） | 低 |

---

## bn_repro.py

实跑：退出码 `1`；stderr 为空。关键 stdout：`Fig4a … → below（论文：尖）`、`Fig4c … → below（论文：宽）`、`与论文一致 → FAIL`、八行命中 `3/8`（Fig4a/4c/4d 命中；Fig4b 判 initial burst、Fig4e 判 accelerating、Fig4f 判 **无发放**、Fig4g/4h 判 initial burst）、`自检总判: 有 FAIL`。
跨目录实跑：退出码 `1`，stderr 末行 `FileNotFoundError: [Errno 2] No such file or directory: '/tmp/bn_adex.py'`。

| 行号 | 类别 | 事实 | 为何是缺陷 | 复现 | 严重度 |
| --- | --- | --- | --- | --- | --- |
| 36-38 | (e) 工作目录 | `adex()` 用 `spec_from_file_location("bn_adex", "bn_adex.py")` 裸文件名。 | 换目录运行即在 `/tmp` 找 `bn_adex.py`，直接 FileNotFoundError；整个复现依赖当前工作目录。 | `cd /tmp && timeout 60 python3 -u /sdcard/dsphworkspace/brainnet-core/bn_repro.py; echo EXIT=$?`（EXIT=1，stderr 报 `/tmp/bn_adex.py` 不存在） | 中 |
| 145-163 | (b) 判据 | 判据 [2] 检验论文自己的两个例子：Fig4a（b=0, V_r=−58）判 below（与论文一致），Fig4c（b=120, V_r=−50）也判 below，而论文要求宽 SAP 的重置点在 V-nullcline 之上。实跑 `与论文一致 → FAIL`，退出码 1。数值：Fig4c `w_reset=121.06 vs w*(V_r)=292.00`。我用真实穿越时刻的 w（`w_cross=1.0605`、`w_reset=121.0605`）复算，仍是 below，排除“取子步起点 w 的插值误差”。 | 复现脚本无法复现论文图注自举的两个例子，判据 [2] 判据/参数不自洽；根因未证实（两种可能：`w*(V_r)` 这个比较面的选取与论文相位图判据不等价，或 bn_adex 的 Fig4c 参数转录有误）。无论哪种，本文件当前输出与论文标签相反。 | `cd /sdcard/dsphworkspace/brainnet-core && timeout 200 python3 -u bn_repro.py; echo EXIT=$?`（`Fig4c … → below（论文：宽）`、`与论文一致 → FAIL`、EXIT=1） | 高 |
| 165-178 | (b) 判据 | 八行逐行复现实跑命中 3/8；`hits >= 5` 的阈值（L178）无来源；`EXPECTED` 把 Fig4b 标 adaptation、Fig4g 标 transient spiking、Fig4h_5 标 irregular spiking，`classify` 分别判成 initial burst / initial burst / initial burst。 | 命中率 3/8 说明分类器的模式判据与论文八行标签大面积不符；阈值 5/8 无依据，且即使达到也只是“过半数”。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_repro.py 2>&1 \| grep "命中"`（命中 3/8） | 中 |
| 118-142 | (c) 度量定义 | 分类阈值为硬编码经验值：`isi[0] < 0.5*mid`（L124）、`x > 3*mid`（L121/L124）、`(isi[-1]-isi[0]) > 0.2*isi[0]`（L131）、`isi[-1] < 0.85*isi[0]`（L135）、`cv > 0.3`（L140）、`n >= 8`（L140）。 | 这些系数没有论文出处，也不是论文给出的判据；论文只对 sharp/broad SAP 与 initial/regular bursting 给了相平面判据，对 delayed/transient/irregular 明确未给解析判据（文件头 L27-29 自认）。用无来源阈值套标签，命中率低是必然。 | 读 L118-142；实跑 `命中 3/8` | 中 |
| 180-186 | (b) 判据 | 混沌段只打印 `轨迹等价=False`，没有 `ok &=`；检验方法是比较两次仿真的脉冲时刻是否逐点相差 `<1e-6`（L184），不是发散率/Lyapunov 量。 | 该段不是可失败判据；“混沌”的论文主张（低维混沌）既未被证实也未被否证，输出措辞（“只报观测”）与文件头 L16 的判据承诺不一致。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_repro.py 2>&1 \| tail -6`（混沌段无 PASS/FAIL） | 中 |
| 67-84 | (a) 运行期 | `reset_side` 自写第二套积分循环（L70-81）来重放轨迹并取放电前的 `ww`（L79），与 `ax.run` 的实现（bn_adex L65-99）不是同一段代码；文件内没有任何断言校验两者给出相同重置点。 | 同一物理轨迹被两套代码分别积分，存在实现漂移风险；本文件没有交叉验证，无法保证 `reset_side` 用的是产生 `spikes` 的那条轨迹。**未证实**是否已实际造成偏差（我用同款积分复算 Fig4c，结论与文件一致）。 | `cd /sdcard/dsphworkspace/brainnet-core && python3 -u bn_repro.py`（reset_side 的判定与 `ax.run` 的脉冲序列无一致性检查） | 低 |

---

## 本组小结

### 缺陷总数与分类计数

共 **47** 条（按“行号/类别/事实/为何缺陷/复现/严重度”表逐行计）。

| 类别 | 条数 | 涉及文件 |
| --- | --- | --- |
| (a) 运行期 | 4 | bn_grow_v2（死代码 1）、bn_grow_time（死代码 1）、bn_perc（死代码 1）、bn_repro（两套积分，未证实 1） |
| (b) 判据缺陷 | 24 | bn_growth 6、bn_grow_time 2、bn_grow_rule 3、bn_spatial 4、bn_perc 2、bn_part 4、bn_repro 3 |
| (c) 度量定义 | 9 | bn_grow_v2 1、bn_grow_time 2、bn_grow_rule 2、bn_spatial 1、bn_perc 2、bn_repro 1 |
| (d) 文件头与代码不一致 | 8 | bn_grow_v2 3、bn_grow_time 1、bn_spatial 1、bn_perc 2、bn_part 1 |
| (e) 依赖当前工作目录 | 2 | bn_grow_v2、bn_repro |

严重度分布：高 5、中 26、低 16。

### 最严重的三条

1. **bn_growth.py L139/L146/L159——正向判据从不进入 `ok`，总判空转。** 实跑对齐比例 0.060 低于同规模随机 0.084、`显著高于随机=False`，文件仍打印 `全部 PASS` 且退出码 0。生长机制在默认参数下（`window=1` vs `T0=8`，延迟恒为 {0,1}）几乎不可能产生 mod 8 对齐，这一实质失败被读出为成功。
2. **bn_part.py L20/L128-143——跨区回路恒为 0，跨区对齐判据所测对象不存在。** 实跑 `回路=42（区内 42 跨区 0）`、`跨区=None`，零分布段直接 `continue`；文件头承诺“区内/跨区各自比较”，实际只兑现一类，仍报全部 PASS。同文件 L13/L18 的“hub→hub 2 拍、至少跨 3 个取值”与实跑 `{1:36, 3:8}`（仅 2 值）矛盾。
3. **bn_grow_time.py L8-14——“全库无显式数值 λ”与论文库事实相反。** `papers/PMC12380633.xml` 明写 `λ_decay = 1/2/10 mm⁻¹`（图 1C）与 `λ_decay = 5 mm⁻¹`（正文）、图 5 真值 1.2–3.2。该错误声明直接导致同一来源下 bn_grow_time 默认 λ=5.117、bn_grow_v2 默认 λ=2.0，两文件默认值冲突。

（另：bn_repro.py 是八件中唯一实跑 `exit=1` 的文件，其 [2] 判据 FAIL 已在 stdout 如实报出，故不列入“最严重三条”的误导性缺陷，但作为唯一失败的复现资产应优先处理。）

### 本组是否有缺陷曾导致错误科学结论

**判断：是（文件输出层已经发生；是否被下游引用进正式结论，未能证实）。**

依据：

- bn_growth.py 已实际输出 `显著高于随机=False` 与 `判据总判: 全部 PASS` 同时出现；若被引用，会得到“候选生长机制能产生对齐回路”的反事实结论（实测 0.060 < 随机 0.084）。
- bn_part.py 已实际输出 `对齐比例：区内=0.0 跨区=None` 与 `判据总判: 全部 PASS`；跨区一半的判据从未执行，“两类比较”的结论不成立。
- bn_perc.py 已把唯一出处 `theory/p45_experiments.py` 的 `f_c = 0.8` 写成 `f_c ≈ 0.85`；头部仍写“给出实测 f_c”，与退役后的实际代码不符。
- bn_grow_time.py 的“λ 数值未知/全库无数值”错误声明，与同库 PMC12380633 的显式 λ 及 bn_grow_v2 的 λ=2.0 有出处声明均冲突。
- 反例（诚实失败，不构成错误结论）：bn_repro.py 自检 FAIL、退出码 1，并在 stdout 明确报出 `与论文一致 → FAIL`、`命中 3/8`；bn_spatial.py 的负例在 `L_MAX=1.8` 下用实测跨区最小距离 2.0406 mm 验证，设置有效。

### 未证实项（不以推测充当结论）

- bn_repro.py L145-163：Fig4c 判 below 的根因未证实。已排除“重置时刻 w 的插值误差”（用真实穿越时刻 `w_cross` 复算仍 below），剩余可能为判据比较面与论文相位图口径不等价、或 `bn_adex.TABLE1` 的 Fig4c 参数转录有误；后者需比对 Naud2008 原表逐列，本轮未完成。
- bn_grow_rule.py：`normalizer` 的积分区间/步长误差对三模型排名的实际影响未做数值敏感性测试。
- bn_spatial.py L145 的 `- 1e-9` 放松在本布局下未触发，是否在其它参数下改变结论未证实。
- 各文件头部引用的“PNAS Nexus pgag285 orchestrator”“R Soc Open Sci 260575 仙女环”“PMC12380633 R=0.93”已在 papers 目录命中对应文本，视为核实；`PMC12319759` 的“skewness=1.71【主控结案 2026-10-09：已在 PMC12319759.xml 原文证实——“The distribution of edge length was highly positively skewed (1.71)”，对象=边长度分布】”本轮只定位到“edge length 高度正偏”的句子，未在该句直接检出 1.71，记为**未证实**。
- 本轮未修改任何被审源码，也未改动 `theory/audit/` 之外的任何文件。

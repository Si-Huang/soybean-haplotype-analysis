---
name: soybean-haplotype-analysis
description: Analyze soybean gene haplotypes with beginner-friendly explanations and consistent scientific figures. Route target genes to promoter analysis and transcription factors to CDS/protein analysis; retrieve SoyOmics data or use audited local inputs, then assess population composition, geography, and available phenotypes.
---

# 大豆基因单倍型分析

把一个大豆基因转化为可复现的分析、图件、源数据和中文解释。面向新手：由助手执行命令、检查环境和解释结果，不要求用户先掌握术语或整理整套输入。

## 入口与分支

- **靶基因 → 启动子**。默认取注释基因边界上游 2 kb，按正负链计算；没有实验证据时称“注释基因边界”，不称已验证 TSS。
- **转录因子 → CDS／蛋白**。确认转录本，CDS 单倍型包括同义与非同义变化；联合重建同一单倍型的编码序列后归并蛋白类型。
- 用户已说明角色就直接执行。角色不明先查可靠注释；“靶基因”是本次调控关系中的身份，不能仅由家族判断。同一 TF 同时作为上游因子的靶基因时按研究问题选择；仍不清楚则用一句通俗问题澄清，不默认两个分支全跑。
- 俗名/Glyma/Wm82 编号先核对版本映射。在线获取器只接受已核实的 **SoyZH13 编号／ZH13.v2**；其他版本用匹配参考版本的本地输入，不改前缀冒充映射。

先说明区域、默认参数和能回答的问题，再执行。默认常见变异 MAF ≥ 0.05，位点缺失率 ≤ 0.10；这是本项目约定，不是普适最佳值。低频功能候选另保留，必要时独立运行 MAF ≥ 0.01 敏感性分析，两次编号互不混用。

## 执行

1. 阅读 [输入与来源](references/inputs-and-sources.md)。定位本 skill 绝对路径；检查 Python 和依赖，缺失时在独立虚拟环境安装 `requirements.txt`。输出放用户项目的新目录。
2. SoyOmics 数据执行 `scripts/fetch_soyomics.py --gene ID --role target|tf --out INPUT_DIR`。可加 `--symbol NAME`、`--transcript ID.m1`、`--gff FILE`。TF 默认先检查 `.m1`，用户指定优先；默认转录本读框冲突时保留审计，按输入指南检查其他完整转录本并说明选择。保存原始响应、参数、时间、哈希、保留位点清单。首次 TF 运行下载官方 GFF，可复用匹配版本本地 GFF。
3. 执行 `scripts/analyze.py --input INPUT_DIR --out RESULT_DIR`。按需传 `--phenotypes CSV`、`--domains CSV`。原始 phased VCF 或已有表格按输入规范转换，不能伪造分相。只有一个类型也如实交付；没有位点且缺少材料分母时报告“当前筛选未保留变异”，不能虚构全体材料类型。
4. 阅读 [统计与扩展](references/analysis.md)。基础包含矩阵、网络、群体频率、映射、报告；可靠坐标才做地图，材料级表型才做表型，候选 TF/PWM 明确才做 motif，环境作为探索模块。缺数据写明缺什么，继续其他独立模块。
5. 遵循 [绘图约定](references/figure-style.md)。核心脚本生成 PDF/SVG/PNG；扩展图复用 `scripts/style.py`。按 [新手讲解](references/beginner-guide.md) 解释结果。
6. 交付前核对编号、分母、材料去重、比例、有效 n 以及图源表一致性。打开图检查拥挤、字号、遮挡和裁切，修正后交付。

所有脚本相对本 skill 根目录，执行时转为绝对路径。README 有可直接复制的练习命令。

## 科学约束

- 区分单倍型内群体组成 `P(population|haplotype)` 与群体内单倍型频率 `P(haplotype|population)`，图上明确分母。
- 材料、二倍体等位拷贝、测定记录是不同单位。拷贝可描述频率，但不作为独立材料产生正式关联 P。内置检验只使用唯一归属的独立材料；杂合另表列出。
- HapSnap 有内部位点过滤、材料缺失/杂合处理、最小拷贝阈值。记录请求/返回清单，不将主要类型当全样本变异谱；不把 `heter=1` 解释为实验分相证明。
- NA 不等于 0、感病或参考型。表型百分比分母为该字段有唯一明确记录的材料 n；同值重复去重，冲突仅排除对应材料×字段。
- 自动蛋白重建限已验证参考 CDS 上明确 SNP；同一密码子 SNP 联合替换。Indel、符号等位、参考不匹配、部分 CDS 时停止自动翻译，保留 CDS 结果并说明需要规范化 VCF/序列重建。称“保留位点定义的蛋白类型”。
- 内置网络为保留位点状态 Hamming 距离最小生成树，边刻度为位点状态差异；不称全序列 TCS/median-joining，不推断祖先或选择方向。
- motif 分数为预测，FIMO 单序列 P/q 不是两等位差异检验。单 SNP 边际频率不能恢复双位点相位。群体频率差异不能直接证明选择、结合或抗病毒因果。
- SMV 公开表型不能替代其他病毒的接种实验。仅依据原始定义解释 R/ER/N/S，ER 不擅称极抗。

## 交付与练习

给出基因身份／角色／版本、完成与缺数据模块、结果及证据边界，链接报告、图和源数据。默认中文报告、英文图文。未运行模块不称完成。

`python scripts/make_demo.py --out DIR` 生成明确标记的合成练习，不代表真实基因。维护脚本时运行 `python -m unittest discover -s tests -v`；该检查覆盖分析不变量和两分支，不替代真实数据审查。

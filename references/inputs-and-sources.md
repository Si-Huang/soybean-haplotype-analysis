# 输入、版本与数据来源

## 最少需要什么

新手先给：基因编号或名字、它在本研究中是靶基因还是 TF、输出位置。已知参考版本、转录本、候选表型、已有材料数据时一起说明。助手先检查项目现有文件，避免让用户重复整理。

在线入口：[SoyOmics](https://ngdc.cncb.ac.cn/soyomics/)。自动获取器只实现 ZH13.v2。俗名或 Glyma 编号不等价于唯一版本；跨版本需数据库对应表或序列／共线性证据，不能字符串换前缀。

获取器保存 gene 查询、链方向、请求位点、最终保留位点、HapSnap 请求和响应、材料信息及 provenance.json。后者包含时间、URL、请求参数、HTTP 状态及 SHA-256。官方 GFF 下载地址在脚本 `GFF_URL`；外部接口不保证长期稳定。`--gff` 本地复用时记录文件哈希，并核对基因、染色体和链方向。

HapSnap 的 count 是返回等位拷贝，可能对缺失、杂合、低频类型继续过滤。必须保留 source_id 和完整请求，不能把只保留的主要类型写成全部种质材料的全变异谱。服务提供的类型和材料分配称为“服务派生”，不宣称实验分相；需要严格相位研究时换用有分相证据的输入。

## 本地输入目录规范

CSV 使用 UTF-8、逗号分隔、首行为列名。ID 保持字符串一致。空值不填 0。

|文件|必需字段|说明|
|---|---|---|
|manifest.json|gene_id, role, assembly, strand, site_columns|role=`target`/`tf`，strand=`+`/`-`；记录区域、阈值、分相来源及转录本|
|sites.csv|site_id, position, alleles|一行一个唯一保留位点，position 为参考基因组 1-based；alleles 如 A/G；启动子另有 relative_position，TF 另有 cds_position|
|haplotypes.csv|haplotype_id, source_id, 每个 site_id 一列|同一行的等位状态必须有共同相位／纯合材料证据；不能把未分相双杂合直接拼接|
|counts.csv|haplotype_id, population, count|每个类型必须完整包含 Soja/Landrace/Improved 三行；无观测写 0，未知不能写 0|
|assignments.csv|sample_id, haplotype_id, population|可选，每行一个等位拷贝；二倍体同型可以两行同 ID，不能超过两行；行数须与 counts 完全一致|
|metadata.csv|sample_id, longitude, latitude|可选，每材料一行；可另有 accession/country/region。sample_id 与输入分配一致，不能直接把展示名当内部 ID|
|reference_cds.fasta|完整单条转录本 CDS|TF 蛋白分支必需，序列为转录本 5′→3′，包含起始密码子、无内部终止；需验证 phase 和坐标|

manifest 的 `site_columns` 顺序必须与 sites.csv 一致。CDS 坐标 `cds_position` 为拼接转录本 1-based，负链基因要反向排列外显子并互补基因组等位。验证后显式填 `cds_phase_valid: true`；否则保留 CDS 分析并跳过蛋白自动重建。

无材料分配时可交付拷贝频率与矩阵，但跳过材料级检验、地理和表型。只有一类群体时必须区分“其他群体确认无观测”与“根本没收集其他群体”；后一情况不得伪造零计数，应先向助手说明并定制分析输入。

## 可选表型与结构域

`--phenotypes CSV` 采用长表：

```csv
sample_id,trait,value
material_1,SMV_original_response,R
material_1,SMV_original_score,3.000
material_2,SMV_original_response,N
```

value 是原始分类或离散等级，`N` 是有效标签而非自动缺失。3、3.0、3.000 合并为同一数值类别；不同 trait 不混合。同一材料×trait 同值重复合并；不同值冲突只排除该字段。连续定量性状、重复环境、重复试验不应直接套类别堆积图，应另按试验设计建模。

表型 sample_id 必须先与材料映射对齐。SoyOmics 的内部 ID 与 PI 等展示编号可能不同，使用 metadata 显式转换并审查一对多映射，不凭行号拼接。随表保存字段字典、物种／病毒／株系、环境、原出处和缺失编码。

`--domains CSV`：`name,start,end,source`。位置为匹配参考转录本蛋白的 1-based 闭区间，source 填可追溯数据库记录／版本／链接。不要把不同转录本或物种的结构域位置直接搬来。

## 原始 VCF 如何接入

本版本未提供通用 VCF→表格转换器。助手可用 [BCFtools](https://samtools.github.io/bcftools/bcftools.html) 等规范化工具，核对参考 FASTA、多等位／Indel 表示、样本和缺失；明确 phased GT、分相块或纯合证据后再建立上述表格。未分相杂合不能仅按 GT 左右顺序构造全基因单倍型。复杂变异的蛋白效应用完整序列重建，不能交给内置 SNP 翻译器。

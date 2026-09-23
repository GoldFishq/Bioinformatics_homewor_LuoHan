# Week 5 Homework 2 — EasyMultiProfiler-Web RNA-seq 完整分析报告

**工具:** EasyMultiProfiler-Web v7（本机部署，前端 http://127.0.0.1:8080 + 后端 API :8000）
**数据:** `tests/RNAseq_output.csv`（19,150 基因 × 24 样本，小鼠转录组）+ `tests/RNAseq_mapping.csv`（6 组 × 4 重复：DMSO / DMSO+LIPUS / T4400 / T4400+LIPUS / T3976 / T3976+LIPUS）
**Session:** `1JTUGYQPrxmsLNy5ivkND6dq`

## 1. 分析流程

| 步骤 | 接口 | 结果 |
|---|---|---|
| 上传数据 | `POST /api/import` | 24 样本 × 19,150 特征，识别为 transcriptomics（assay=counts） |
| 主比较差异分析 | `POST /api/workflows/transcriptomics/analyze/differential`（method=DESeq2, Group: T4400 vs DMSO, filter_low=TRUE, 每组 n=4） | 16,760 基因进入检验 |
| PCA（样本 QC） | `POST /api/analyze/dimension`（method=PCA） | 两组沿 PC2 分离，无离群样本 |
| 火山图 | `POST /api/workflows/transcriptomics/visualize/volcano` | `EMP_volcano_T4400_vs_DMSO.png` |
| 热图 | `POST /api/workflows/transcriptomics/visualize/heatmap` | `EMP_heatmap_T4400_vs_DMSO.png`（top-66 高变特征，z-score） |
| 结果导出 | 本文件同目录 CSV | `EMP_DESeq2_T4400_vs_DMSO_results.csv`（16,760 行全表） |

## 2. 主要结果（DESeq2，padj < 0.05 且 |log2FC| ≥ 1）

- **上调：约 168–172 个基因**（全表严格 fdr<0.05 计 168 个；火山图按 padj≤0.05 计 172 个）
- **下调：约 66–68 个基因**
- 方向已用原始 count 独立验证：Cep290（DMSO 均值 45.2 → T4400 124.8）、Il1a（21.8 → 124.5）、Dcn（2043 → 6225），均为 T4400 组更高，即**正 log2FC = 在 T4400 中上调**。

### 上调基因中的亮点（见火山图标注）
- **炎症/细胞因子通路：** Il1a、Il6、Ccl3、Ptgs2（COX-2）、Cxcl10 —— 先天免疫与炎症反应激活
- **ECM/成纤维细胞活化：** Col1a1、Col1a2、Col6a1、Dcn、Fn1、Sparc、Thbs4 —— 基质重塑与纤维化相关信号
- **ER 应激/其他：** Ero1l、Axl、Havcr2（KIM-1，损伤标志物）、Itgb7

### 下调基因中的亮点
- Ptprv、Mest、Cdkn1c（印记基因，细胞增殖调控）、Stra6（视黄醇转运）、Lrrc15、Ucma

## 3. 生物学解读

在 DMSO 对照背景下，化合物 T4400 处理诱导了一个以**炎症应答 + 细胞外基质/成纤维细胞活化**为核心的转录程序：Il1a、Il6、Ccl3 等细胞因子和 Cxcl10 趋势显著上调，同时大量胶原基因（Col1a1/Col1a2/Col6a1）与 Dcn、Fn1、Sparc 等基质基因协同上调（热图中 ECM 模块清晰分层）。这一模式提示 T4400 可能触发了组织损伤–修复反应或纤维化样应答，Havcr2（KIM-1）作为肾小管损伤经典标志物的上调进一步支持存在组织应激。下调基因以增殖/印记调控（Cdkn1c）与代谢转运（Stra6）为主，与上调的损伤应答程序互补。

**局限性：** ① 每组仅 n=4，检验效力有限——DMSO+LIPUS vs DMSO 比较中最小 p 值仅 2e-4，BH 校正后无基因存活（fdr 全为 1 是统计正确的，不是 bug）；② 富集分析依赖的 Reactome 数据库包在本机无法安装（无 Rtools），KEGG/GO 富集未在本次报告中展示；③ log2FC 未经 apeglm 收缩，效应量估计在小样本下偏乐观；④ 结论需 Western blot/qPCR 等实验验证。

## 4. 与 Homework 1 概念的衔接

- 同样使用 **DESeq2 + 原始整数 count 矩阵**（未做 TPM/CPM 预转换）；
- 同样遵循**先过滤低表达基因**（filter_low）→ 检验 → **padj 阈值 + 效应量阈值联合判定**的流程；
- Homework 1 中额外做了 apeglm LFC 收缩与 `~ batch + condition` 设计（该数据无批次混杂），Homework 2 由 EMP-Web 平台一键完成 DESeq2 检验，并通过 API 驱动全流程，二者互为印证。

## 5. 排障记录（过程中修复的环境问题）

1. 后端启动失败：R 缺 `plumber`、`metap`（其依赖 `multtest` 缺失）→ 已安装；
2. **tidybulk 1.19.2（devel 版）存在 bug**：`.scale_abundance_se` 引用未定义变量 `scaled_string`，导致所有差异分析方法报 `object 'scaled_string' not found` → 升级到 release 版 tidybulk 2.0.1 后解决；
3. 热图依赖 `pheatmap` → 已安装。

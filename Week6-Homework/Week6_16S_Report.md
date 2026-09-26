# Week 6 — 16S Microbiome Analysis Report

**Course:** Bioinformatics SUAT 2026 Fall (EMP2026)
**Student:** Luo Han (SUAT24000110) · GitHub: GoldFishq / Bioinformatics_homewor_LuoHan
**Tool:** EasyMultiProfiler-Web (EMP-web) v9.0.4 — backend engine run in-process
**Date:** 2026-09-23
**Assignment:** Use the 16S data in the EMP `tests/` folder to run the full microbiome
analysis on EMP-web, then generate a data-driven scientific hypothesis and submit
via the EMP-web system.

---

## 1. Data & Design

| Item | Value |
|---|---|
| Input taxonomy | `tests/16S_level-7.csv` (species-level, `k__;p__;...;s__` sep `;`) |
| Metadata | `tests/16S_mapping.csv` (SampleID, Group, Group_sub) |
| Raw features | 470 species-level taxa × 132 sample columns |
| Grouped samples | **130** (2 `K_XYL_*` columns have no Group → dropped) |
| Patients | **66**, longitudinal: each has a `_01` (before) and `_02` (after) sample |
| Groups | IBS_before 36, IBS_after 36, UC_before 29, UC_after 29 |
| Paired sets | IBS 35 pairs, UC 29 pairs |
| Collapsed to | **137 genera** (Genus level, all kept) |

The design is **paired longitudinal** (same patient before vs after treatment), so
all inferential tests below use **paired** statistics matched by patient.

---

## 2. Workflow (EMP-web 16S pipeline)

`import (tax, Species, ";")` → `prepare taxonomy → Genus` →
`alpha diversity` → `beta diversity (Bray–Curtis PCoA/PCA/NMDS)` →
`visualize (boxplot / scatter / heatmap / barplot)` → `differential taxa (wilcox)` →
`interpret`. The same backend helpers the EMP-web Plumber API loads were executed
in a standalone R session (the live API has a JSON-persistence bug on import; the
engine itself is identical). The "interpret" statistics (paired Wilcoxon, PERMANOVA,
effect sizes) were computed on the prepared Genus matrix.

---

## 3. Results

### 3.1 Alpha diversity (Shannon, paired Wilcoxon signed-rank)

| Disease | before (mean) | after (mean) | median Δ | p |
|---|---|---|---|---|
| IBS | 1.664 | 1.589 | −0.079 | **0.177** (ns) |
| UC  | 1.616 | 1.468 | −0.121 | **0.023** (↓ significant) |

- Observed species: IBS p=0.62, UC p=0.98 (both ns) → the Shannon drop is an
  **evenness** effect, not a loss of species count.
- Simpson: IBS p=0.36, UC p=0.092 (trend, same direction).

### 3.2 Beta diversity (Bray–Curtis)

| Model | R² | F | p |
|---|---|---|---|
| Group | 0.032 | 1.38 | 0.149 |
| Disease (IBS vs UC) | 0.024 | 3.12 | **0.008** |
| Time (before vs after) | 0.005 | 0.69 | **0.671** |
| Disease×Time | 0.032 | 1.38 | 0.144 |

- **Dispersion homogeneity (betadisper):** F=0.157, p=0.93 → groups have equal
  spreads, so the PERMANOVA p-values are valid (no dispersion confound).
- **Within-patient Bray–Curtis = 0.345**, i.e. **0.75× the population mean** → each
  patient's microbiota is more similar to itself over time than to other patients.
- Interpretation: treatment does **not** restructure the community at genus resolution;
  IBS and UC communities are intrinsically different at baseline (Disease p=0.008).

### 3.3 Differential abundant genera (paired Wilcoxon on relative abundance, BH-corrected)

- **BH q < 0.05 genera: 0** in both IBS and UC → no single genus is a decisive
  treatment marker at genus resolution.
- Nominal candidates (p < 0.05, exploratory only):
  - **IBS (mostly gain):** *A. muciniphila* +1.92 (p=0.013), *L. reuteri* +4.27
    (p=0.009), *B. animalis* +5.22 (p=0.009); trending down: *B. uniformis*,
    *B. ovatus*, *B. gnavus*.
  - **UC (mostly loss):** *P. parvula* −2.35 (p=0.019), *B. caccae* −2.81
    (p=0.024), *B. mucilaginosa* −1.19 (p=0.025), *B. fragilis* −1.73 (p=0.046),
    *B. uniformis* −0.86 (p=0.05).

---

## 4. Scientific Hypothesis

> **H:** In ulcerative colitis (UC) — but not in IBS — the therapeutic intervention
> reduces gut-microbiome **alpha diversity (Shannon)** without detectably reshaping
> overall community structure at genus resolution. The treatment acts as an
> *evenness eroder* rather than a *community replacer*: diversity falls because
> several moderate-abundance taxa are eroded in a distributed way, not because a
> dominant keystone genus is replaced.

**Parameters used to support the hypothesis**

1. **Alpha diversity — paired Wilcoxon signed-rank on Shannon index:**
   UC before→after p=0.023, Δmedian=−0.12 (significant decrease); IBS p=0.18 (ns).
   Observed-species paired p=0.98 (UC) shows the change is in *evenness*, not
   species richness → consistent with loss of rare/low-abundance taxa.
2. **Beta diversity — Bray–Curtis PCoA + PERMANOVA:** Time term R²=0.005 (p=0.67,
   ns) and homogeneous dispersion (p=0.93) confirm **no community-level restructuring**;
   within-patient Bray (0.345) = 0.75× population spread → patients stay stable.
3. **Differential abundance — paired Wilcoxon on genus relative abundance, BH-corrected:**
   0 FDR-significant genera in either disease, with only distributed nominal candidates
   (UC: *P. parvula, B. caccae, B. fragilis, B. uniformis* trending down; IBS:
   *A. muciniphila, L. reuteri* trending up) → no keystone replacement, supporting a
   distributed evenness erosion.
4. **Stratification check — Disease PERMANOVA p=0.008:** IBS and UC communities differ
   at baseline, justifying disease-stratified (not pooled) analysis.

**Predicted follow-up:** if the hypothesis is correct, (a) species/ASV- or
metagenomic-functional resolution will localize the eroded taxa to butyrate producers
and mucus-associated commensals; (b) the UC diversity loss will correlate with
clinical response (Group_sub poor/great), testable with a larger paired cohort.

---

## 5. Limitations

- **Compositional data:** relative abundances; Bray–Curtis/PERMANOVA on proportions can
  be biased by compositionality. Conclusions are descriptive at genus resolution.
- **Multiple testing:** no taxon survives FDR → candidate genera are exploratory.
- **Power:** 35 (IBS) / 29 (UC) pairs; subtle effects may be underpowered.
- **Confounders:** `Group_sub` (poor/great) showed no baseline alpha difference, but
  metadata lacks a site/batch variable ("Site completely confounded with treatment" is
  a generic caution from the reading material and cannot be fully ruled out here).
- **Resolution:** genus-level; species/ASV or functional profiling may reveal mechanism.

---

## 6. Reproducibility

- `analysis/run_emp_16s_inproc.R` — EMP-web import + one-click `run_all_m16s` bundle.
- `analysis/analyze_week6_paired.R` — paired alpha / PERMANOVA / differential statistics.
- Inputs copied to `EMP2026/Week_06/microbiome_16s/weekly/runs/<run_id>/data/`.
- Result tables + plots in `.../results/` and `.../plots/`; this report in `.../teaching/`.

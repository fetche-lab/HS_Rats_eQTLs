# eQTL QC, Candidate Prioritization & Trans-Band Analysis

**Cohort:** HSNIH-Palmer outbred rat population (GeneNetwork)  
**Tissues:** Adipose & Liver RNA-Seq  
**Author:** Felix Lisso  

This repository contains the complete analysis pipeline, results, and reporting for a quality-control and prioritization study of expression quantitative trait loci (eQTLs) in rat adipose and liver tissues. The project intentionally diverges from the standard cis-only approach (Hong-Le et al. 2023) by retaining and systematically investigating trans-eQTL signals.

---

## Repository Structure

```
.
├── data/           # Raw input data
├── scripts/        # Python analysis pipeline
├── results/        # Generated CSV tables
├── figures/        # Publication-ready plots
├── reports/        # LaTeX reports, presentations & markdown summaries
├── docs/           # Instructions and plot guides
└── hs_tmp/         # Legacy / intermediate / temp files (not part of main deliverables)
```

---

## Quick Start

1. **Data** is in `data/`.
2. **Run the pipeline** in order:
   ```bash
   python scripts/eqtl_analysis_pipeline.py
   python scripts/cross_tissue_analysis.py
   python scripts/generate_final_candidates.py
   ```
3. **View outputs** in `results/` and `figures/`.
4. **Read the full report** in `reports/progress_report01.pdf`.
5. **View the 10-minute presentation** in `reports/presentation.pdf`.

---

## Folder Details

### `data/` — Raw Input Data

| File | Description |
|------|-------------|
| `adipose_rlog_table.csv` | GeneNetwork rlog-transformed adipose RNA-Seq expression matrix (HSNIH-Palmer r4) |
| `liver_rlog_table.csv` | GeneNetwork rlog-transformed liver RNA-Seq expression matrix (HSNIH-Palmer r4) |

These files are pre-filtered eQTL summary tables from [GeneNetwork 2](https://genenetwork.org). Each row represents the strongest association peak for one gene, with columns for chromosome, position, `-log10(P)`, and effect size.

---

### `scripts/` — Analysis Pipeline

| File | Purpose |
|------|---------|
| `eqtl_analysis_pipeline.py` | **Main pipeline.** Loads raw data, classifies cis vs trans (4 Mb window), applies four trans flags (`MEGA_EFFECT`, `PSEUDOGENE_TRANS`, `LOC_UNCHARACTERIZED`, `CROSS_TISSUE_TRANS`), assigns Tier 1/2/3 scores, and outputs annotated tissue-specific CSVs. |
| `cross_tissue_analysis.py` | **Cross-tissue integration.** Merges adipose and liver results, computes master scores, detects 5 Mb trans-bands, identifies 33 cross-tissue trans hotspots, and generates all cross-tissue comparison plots. |
| `generate_final_candidates.py` | **Candidate prioritization.** Produces Tier 1 lists, pangenomic candidate lists (529 flagged/uncharacterized genes), and the final actionable candidate table. |
| `analyze_eqtls.py` | Early exploratory script for basic cis/trans classification and summary statistics. Superseded by `eqtl_analysis_pipeline.py` but retained for provenance. |

**Execution order:** `eqtl_analysis_pipeline.py` → `cross_tissue_analysis.py` → `generate_final_candidates.py`

---

### `results/` — Output Tables

| File | Description |
|------|-------------|
| `Adipose_annotated_eqtls.csv` | Full adipose eQTL table with `cis_trans`, `trans_flag`, `tier`, and tissue-specific columns (891 genes). |
| `Liver_annotated_eqtls.csv` | Full liver eQTL table with annotations (570 genes). |
| `Adipose_Tier1_candidates.csv` | 118 high-confidence adipose cis-eQTLs (Tier 1). |
| `Liver_Tier1_candidates.csv` | 116 high-confidence liver cis-eQTLs (Tier 1). |
| `cross_tissue_merged.csv` | 1,264 merged genes with `pleiotropy_class` and `master_score`. |
| `cross_tissue_master_candidates.csv` | Top-scoring candidates from cross-tissue integration. |
| `trans_bands.csv` | **73 trans-bands** identified via 5 Mb binning, with `N_Genes`, `Max_logP`, and `Flags`. |
| `pangenomic_candidates.csv` | 529 flagged/uncharacterized genes (LOC/pseudogene/trans) proposed for pangenomic re-mapping. |
| `final_actionable_candidates.csv` | Final prioritized list combining Tier 1 cis genes and trans-band hotspots for follow-up. |

**Key columns in annotated files:**
- `cis_trans`: Classification based on 4 Mb window
- `trans_flag`: Granular flag string (`MEGA_EFFECT`, `PSEUDOGENE_TRANS`, `LOC_UNCHARACTERIZED`, `CROSS_TISSUE_TRANS`)
- `tier`: `Tier_1_HighConfidence`, `Tier_2_Validate`, or `Tier_3_FlaggedInspect`

---

### `figures/` — Publication-Ready Visualizations

#### Quality Control & Overview
| File | What it shows |
|------|---------------|
| `eQTL_QC_summary.png` | **Four-panel QC figure:** (1) cis vs trans bar chart, (2) tier distribution, (3) effect size violin by class, (4) trans-flag breakdown. |
| `eQTL_manhattan.png` | Genome-wide Manhattan plot of adipose and liver eQTLs, with cis peaks in orange/blue and trans peaks in red. |

#### Diagnostic Plots (Pre-filtered Data)
| File | What it shows |
|------|---------------|
| `ecdf_cis_trans.png` | **Survival ECDF** of `-log10(P)` for cis vs trans hits. Validates that trans hits have a fatter tail of extreme p-values despite being pre-filtered. |
| `violin_cis_trans.png` | Violin-box-strip plot comparing effect-size and significance distributions between cis and trans classes. |

#### Cross-Tissue Comparisons
| File | What it shows |
|------|---------------|
| `cross_tissue_effectsize_comparison.png` | Adipose vs liver effect sizes for 197 shared genes. Black triangles mark **33 cross-tissue trans hotspots**. |
| `cross_tissue_logp_comparison.png` | Adipose vs liver `-log10(P)` for shared genes. |
| `cross_tissue_pleiotropy.png` | Pie chart of pleiotropy classes: conserved cis (135), flipped cis (28), trans both (33), tissue-specific (1). |
| `cross_tissue_specificity.png` | Tissue specificity scores for top master candidates. |
| `cross_tissue_genomic_context.png` | Genomic distribution of cross-tissue hits by chromosome. |
| `cross_tissue_tier1_overlap.png` | Venn diagram of Tier 1 adipose vs liver overlap. |

#### Trans-Band Analysis
| File | What it shows |
|------|---------------|
| `trans_band_lollipop.png` | **Lollipop plot** ranking the top 20 trans-bands by gene count and max `-log10(P)`. |
| `trans_band_landscape.png` | **Genome-wide landscape** of all trans-bands across chromosomes. Band width = gene count; height = max significance. |
| `trans_band_faceted.png` | Per-chromosome faceted panels showing trans-band density, eliminating x-axis overlap. |
| `chr7_zoom.png` | **Zoomed-in Chr7:0–5 Mb** trans-band showing 9 uncharacterized LOC genes, with the 0–5 Mb window highlighted. |

---

### `reports/` — Written Outputs

| File | Description |
|------|-------------|
| `progress_report01.pdf` | **13-page formal progress report** (LaTeX). Covers pipeline, key findings, divergence from Hong-Le et al. 2023, and next steps. |
| `progress_report01.tex` | Source for the progress report. |
| `presentation.pdf` | **10-slide Beamer presentation** for group meeting (16:9 widescreen). |
| `presentation.tex` | Source for the presentation. |
| `eQTL_analysis_report.md` | Markdown summary of tissue-specific eQTL analysis (auto-generated by pipeline). |
| `cross_tissue_report.md` | Markdown summary of cross-tissue integration and trans-band results. |
| `final_recommendations.md` | Actionable recommendations: Tier 1 genes to validate, trans-bands to inspect, pangenomic proposal. |
| `supervisor_response.md` | Detailed response to supervisor feedback, including trans-band table and flagging rationale. |

---

### `docs/` — Project Notes

| File | Description |
|------|-------------|
| `instructions_pj.org` | Original project instructions and requirements. |
| `plot_guide.md` | Guidelines for figure generation and visual standards. |

---

### `hs_tmp/` — Legacy & Temporary Files

This folder holds older drafts, intermediate bundles, LaTeX build artifacts, and superseded outputs. **Not required for reproducibility.**

Examples:
- `db220252.pdf` — Hong-Le et al. 2023 literature PDF
- `progress_report.tex/.pdf` — Earlier draft of the progress report
- `progress_report01_bundle.zip` — Old zip archive of bundled deliverables
- `final_edits/` — Complete copy of files from an intermediate bundling step
- `reporting/` — Broken / partial PDF builds from an earlier compilation attempt
- `qqplot_cis_trans.png` — Superseded by `ecdf_cis_trans.png` and `violin_cis_trans.png`
- `.aux`, `.log`, `.out`, `.nav`, `.snm`, `.toc` — LaTeX auxiliary files

---

## Key Findings (Summary)

1. **Retain trans-eQTLs, don't discard.** All 157 trans hits are kept with four granular annotation flags.
2. **33 cross-tissue trans hotspots.** These replicate in both adipose and liver, arguing against random noise.
3. **7 trans-bands identified** via 5 Mb clustering; 3 replicate across tissues (Chr8, Chr5, Chr4).
4. **Chr7:0–5 Mb is the strongest band** — 9 uncharacterized LOC genes clustered at the telomere.
5. **135/197 shared genes** have conserved cis regulation across tissues (same peak, same direction).
6. **Intentional divergence** from Hong-Le et al. 2023: 4 Mb cis window, extreme-tail focus, trans-signal mining.

---

## Methods at a Glance

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Cis window | 4 Mb | Captures distal enhancers within topologically associating domains (TADs) |
| Trans-band bin | 5 Mb | Balances resolution with statistical power for clustering |
| Pre-filter threshold | `-log10(P) ≥ 20` | Input data from GeneNetwork extreme tail |
| Cross-tissue overlap | Symbol-based | Per-gene best peak merged across tissues |
| Master score | Composite | `-logP` + effect-size consistency + cross-tissue replication |

---

## Next Steps (From Final Recommendations)

**Track 1 — Cis Candidates**
- Verify symbols and check peaks in RGD / Ensembl
- Colocalization with metabolic GWAS (COLOC)
- Fine-mapping on top 5–10 loci

**Track 2 — Trans-Bands**
- Raw read-depth analysis at top 3 bands
- Founder haplotype tracing (8 HS founders)
- Pangenomic re-mapping (Minigraph-Cactus with founder assemblies)

**Suggestive Hypothesis:** Population-segregating structural variants, absent from the linear reference, act as master regulators of multiple expression traits.

---

## Dependencies

- Python ≥ 3.9
- pandas
- numpy
- matplotlib
- seaborn
- scipy

LaTeX compilation requires a standard TeX Live distribution with `beamer`, `booktabs`, `natbib`, and `xcolor`.

---

## Citation / Reference

> Hong-Le et al. (2023). *Expression Quantitative Trait Loci in Adipose and Liver Tissues From an HS Rattus norvegicus Population.* Diabetes, 72:135–148.  
> `db220252.pdf` (available in `hs_tmp/`)

# HSNIH-Palmer Adipose / Liver eQTL — QC, Candidate Prioritization & Trans-Band Analysis

Quality control, tiered scoring, cross-tissue integration, and trans-band clustering for
expression quantitative trait loci mapped in adipose and liver from the **HSNIH-Palmer**
outbred rat population. The package intentionally retains and investigates trans-eQTL signals
rather than discarding them as artifacts, diverging from the cis-only convention in the
published literature.

This README documents the full chain of custody — **who, when, where, what** — for both the
upstream pre-filtered eQTL tables (from GeneNetwork) and the downstream analysis in this
package, so that every number on every figure can be traced back to a source.

---

## 1. The data

| Field | Value |
|---|---|
| **Cohort** | HSNIH-Palmer outbred rats (*Rattus norvegicus*) |
| **Tissues** | Adipose (subcutaneous / visceral) and Liver |
| **Platform** | RNA-Seq, rlog-transformed |
| **Pre-filter** | **−log₁₀(P) ≥ 20** (extreme tail only) |
| **Adipose genes** | **891** (790 cis, 101 trans) |
| **Liver genes** | **570** (514 cis, 56 trans) |
| **Reference genome (GN2)** | **mRatBN7.2** |
| **Reference genome (Hong-Le 2023)** | Rnor_6.0 |
| **Trait page** | https://genenetwork.org (search: *HSNIH-Palmer Adipose RNA-Seq*, *HSNIH-Palmer Liver RNA-Seq*) |

> **Note on the input and scope.** The starting material is not raw reads or counts, but
> **GeneNetwork pre-filtered eQTL summary tables** — one row per gene, reporting only the
> strongest association peak (chromosome, position, −logP, effect size). The genome-wide
> background below −logP = 20 has already been truncated by the upstream pipeline.
>
> Consequently, this package focuses on **evaluating and re-analyzing the downloaded trait
> values** (peak assignments, cis/trans classification, cross-tissue integration, trans-band
> clustering) rather than re-mapping from raw data or regenerating QTL plots from scratch.
> For guidance on viewing the original QTL plots and interactive trait pages directly in
> GeneNetwork 2, see `docs/gn2_trait_visualization_guide.md`.

---

## 2. Provenance of the upstream data (GeneNetwork)

All of the following was read from the GeneNetwork 2 download interface and the file metadata
embedded in the rlog tables; nothing here is inferred unless explicitly marked.

**WHO** — mapping and pre-filtering performed by the GeneNetwork 2 compute pipeline;
data set assembled by **Prof Pjotr Prins** / GeneNetwork team. Run initiated by
**Felix Lisso**.

**WHEN** — rlog tables downloaded **2026-05-26** by Felix Lisso. Upstream GEMMA mapping
run timestamp **[TO BE CONFIRMED from GeneNetwork metadata]**.

**WHERE** — GeneNetwork 2 production server (`genenetwork.org`). Compute host details
not independently documented in the downloaded tables.

**WHAT** — RNA-Seq expression QTL mapping for HSNIH-Palmer rats. Because the raw
associations are not included, the exact tool version, covariates, kinship model, and
permutation parameters used by GeneNetwork are **not independently verifiable from these
files**.

| Parameter | Value | Source |
|---|---|---|
| GeneNetwork dataset (adipose) | `HSNIH-Palmer Adipose RNA-Seq (Feb26) rlog` | GN2 dataset name |
| GeneNetwork dataset (liver) | `HSNIH-Palmer Liver RNA-Seq (Feb26) rlog` | GN2 dataset name |
| Genotypes | `HSNIH-Palmer_r4 Genotypes` | GN2 dataset name |
| Pre-filter threshold | −log₁₀(P) ≥ 20 | file content; max −logP = 169.4 (adipose), 175.6 (liver) |
| Genotype build | **HSNIH-Palmer_r4** | file name |
| Reference genome | **mRatBN7.2** | GeneNetwork 2 metadata |
| Mapping tool | **GEMMA** | GeneNetwork 2 standard pipeline |

> **Two upstream metadata gaps — flagged, not hidden:**
> 1. The **exact GEMMA version** (assumed v0.98.5 based on GeneNetwork production at the
>    time), the **kinship model** (leave-one-chromosome-out or full GRM), and the exact
>    **covariates** are not recorded in the downloaded tables.
> 2. The **genome-wide significance threshold** used by GeneNetwork to pre-filter to
>    −logP ≥ 20 is not stated; it may be a hard cut-off, an empirical permutation threshold,
>    or a Bayes-factor heuristic. We treat −logP = 20 as a given floor, not an FDR-corrected
>    genome-wide threshold.

---

## 3. Provenance of this analysis package (DOWNSTREAM)

**WHO** — **Felix Lisso**, supervised by **Prof Pjotr Prins**, GeneNetwork Lab, UTHSC.

**WHEN** — downstream analysis performed **2026-05-28 to 2026-06-03**.

**WHERE** — local Linux workstation. Software: Python 3, `pandas` 2.x, `numpy`, `matplotlib`
3.x, `seaborn`, `scipy`. LaTeX compilation via TeX Live 2023 (`pdflatex`).

**WHAT** — for the pre-filtered adipose and liver eQTL tables:

1. **Cis vs trans classification** with a **4 Mb window** (not the 1 Mb used by Hong-Le
   et al. 2023). Genes whose peak lies within 4 Mb of their own transcription start are
   called cis; all others are trans.
2. **Granular trans flags** — replaced the previous boolean `flag_artifact` with four
   categorical flags: `MEGA_EFFECT`, `PSEUDOGENE_TRANS`, `LOC_UNCHARACTERIZED`,
   `CROSS_TISSUE_TRANS`. All 157 trans hits are retained.
3. **Three-tier scoring** — `Tier_1_HighConfidence`, `Tier_2_Validate`,
   `Tier_3_FlaggedInspect`.
4. **Cross-tissue integration** — merged 1,264 genes by symbol, computed master scores,
   and assigned pleiotropy classes (`shared_conserved_cis`, `shared_cis_flipped`,
   `trans_both_hotspot`, `tissue_specific_cis`).
5. **Trans-band clustering** — binned trans peaks into 5 Mb windows, yielding **73 bands**;
   **7 bands** exceeded a heuristic enrichment threshold.
6. **Visualization** — QC summary, Manhattan plot, cross-tissue effect-size / −logP
   comparisons, pleiotropy pie chart, trans-band lollipop / landscape / faceted plots,
   Chr7 zoom, ECDF and violin diagnostics.
7. **Reporting** — 13-page progress report (`progress_report01.pdf`) and 10-slide Beamer
   presentation (`presentation.pdf`).

---

## 4. Headline results

- **162 Tier 1 high-confidence cis-eQTLs** (118 adipose + 116 liver, 72 overlapping).
- **157 trans-eQTLs retained**, not discarded — 33 replicate across both tissues
  (cross-tissue trans hotspots, marked with black triangles on effect-size plots).
- **7 trans-bands** identified by 5 Mb clustering; **3 replicate across tissues**
  (Chr8:40–45 Mb, Chr5:160–165 Mb, Chr4:105–115 Mb).
- **Chr7:0–5 Mb** is the strongest band — **9 uncharacterized LOC genes** clustered at the
  telomere, proposed for read-depth / CNV inspection.
- **135 / 197 shared genes** have conserved cis regulation (same peak, same direction) across
  adipose and liver.
- ~**40 % of all hits** are uncharacterized LOC symbols or pseudogenes — retained as
  pangenomic candidates rather than filtered out.

**This is an intentional divergence from Hong-Le et al. 2023**, who used a 1 Mb cis window,
reported genome-wide at FDR < 0.001, and essentially ignored trans-eQTLs. We focus on the
extreme tail, use a 4 Mb cis window, and systematically mine the trans signal.

---

## 5. File manifest

### Data
| File | What it is |
|---|---|
| `data/adipose_rlog_table.csv` | GeneNetwork pre-filtered adipose eQTL summary table. Source: **HSNIH-Palmer Adipose RNA-Seq (Feb26) rlog**. |
| `data/liver_rlog_table.csv` | GeneNetwork pre-filtered liver eQTL summary table. Source: **HSNIH-Palmer Liver RNA-Seq (Feb26) rlog**. |

### Scripts
| File | What it is |
|---|---|
| `scripts/eqtl_analysis_pipeline.py` | Main pipeline. Loads raw data, cis/trans classification (4 Mb), trans-flag assignment, tier scoring. Outputs annotated tissue-specific CSVs. |
| `scripts/cross_tissue_analysis.py` | Cross-tissue integration. Merges adipose + liver, computes master scores, detects 5 Mb trans-bands, identifies 33 cross-tissue trans hotspots. Generates all cross-tissue comparison plots. |
| `scripts/generate_final_candidates.py` | Candidate prioritization. Produces Tier 1 lists, pangenomic candidate list (529 flagged/uncharacterized genes), and final actionable candidate table. |
| `scripts/analyze_eqtls.py` | Early exploratory script for basic cis/trans classification. Superseded by `eqtl_analysis_pipeline.py` but retained for provenance. |

### Results
| File | What it is |
|---|---|
| `results/Adipose_annotated_eqtls.csv` | Full adipose eQTL table (891 genes) with `cis_trans`, `trans_flag`, `tier`. |
| `results/Liver_annotated_eqtls.csv` | Full liver eQTL table (570 genes) with annotations. |
| `results/Adipose_Tier1_candidates.csv` | 118 high-confidence adipose cis-eQTLs. |
| `results/Liver_Tier1_candidates.csv` | 116 high-confidence liver cis-eQTLs. |
| `results/cross_tissue_merged.csv` | 1,264 merged genes with `pleiotropy_class` and `master_score`. |
| `results/cross_tissue_master_candidates.csv` | Top-scoring candidates from cross-tissue integration. |
| `results/trans_bands.csv` | 73 trans-bands from 5 Mb binning, with `N_Genes`, `Max_logP`, `Flags`. |
| `results/pangenomic_candidates.csv` | 529 flagged / uncharacterized genes proposed for pangenomic re-mapping. |
| `results/final_actionable_candidates.csv` | Final prioritized list combining Tier 1 cis genes and trans-band hotspots. |

### Figures
| File | What it is |
|---|---|
| `figures/eQTL_QC_summary.png` | Four-panel QC: cis vs trans counts, tier distribution, effect-size violin, trans-flag breakdown. |
| `figures/eQTL_manhattan.png` | Genome-wide Manhattan plot (adipose + liver). |
| `figures/ecdf_cis_trans.png` | Survival ECDF of −logP for cis vs trans hits. |
| `figures/violin_cis_trans.png` | Violin-box-strip plot of effect size and significance by class. |
| `figures/cross_tissue_effectsize_comparison.png` | Adipose vs liver effect sizes; black triangles = 33 cross-tissue trans hotspots. |
| `figures/cross_tissue_logp_comparison.png` | Adipose vs liver −logP for shared genes. |
| `figures/cross_tissue_pleiotropy.png` | Pleiotropy class pie chart (135 conserved cis, 28 flipped, 33 trans both, 1 tissue-specific). |
| `figures/cross_tissue_specificity.png` | Tissue specificity scores for top master candidates. |
| `figures/cross_tissue_genomic_context.png` | Genomic distribution of cross-tissue hits by chromosome. |
| `figures/cross_tissue_tier1_overlap.png` | Venn diagram of Tier 1 adipose vs liver overlap. |
| `figures/trans_band_lollipop.png` | Lollipop plot ranking top 20 trans-bands by gene count and max −logP. |
| `figures/trans_band_landscape.png` | Genome-wide trans-band landscape (width = gene count, height = significance). |
| `figures/trans_band_faceted.png` | Per-chromosome faceted panels of trans-band density. |
| `figures/chr7_zoom.png` | Zoomed Chr7:0–5 Mb trans-band with 9 uncharacterized LOC genes highlighted. |

### Reports
| File | What it is |
|---|---|
| `reports/progress_report01.pdf` | 13-page formal progress report (LaTeX). |
| `reports/progress_report01.tex` | Source for the progress report. |
| `reports/presentation.pdf` | 10-slide Beamer presentation (16:9 widescreen). |
| `reports/presentation.tex` | Source for the presentation. |
| `reports/eQTL_analysis_report.md` | Auto-generated markdown summary of tissue-specific eQTL analysis. |
| `reports/cross_tissue_report.md` | Auto-generated markdown summary of cross-tissue integration. |
| `reports/final_recommendations.md` | Actionable recommendations: Tier 1 validation, trans-band inspection, pangenomic proposal. |
| `reports/supervisor_response.md` | Detailed response to supervisor feedback, including trans-band table and flagging rationale. |

### Documentation
| File | What it is |
|---|---|
| `docs/instructions_pj.org` | Original project instructions and requirements. |
| `docs/plot_guide.md` | Visual standards and figure-generation guidelines. |
| `docs/gn2_trait_visualization_guide.md` | How to view original QTL plots, founder haplotypes, and allele-effect diagrams for candidate genes directly in **GeneNetwork 2**. Complements the downloaded-table analysis by linking results back to the interactive source. |

---

## 6. Known limits (so they cannot be used against us)

- **Pre-filtered input.** The genome-wide background below −logP = 20 is **already missing**.
  We cannot compute a conventional genomic inflation factor (λ), permutation thresholds, or
  empirical FDRs from these tables alone. The ECDF and violin plots are therefore
  descriptive, not inferential.
- **4 Mb cis window is a deliberate divergence.** The literature (Hong-Le et al. 2023) used
  1 Mb. Our 4 Mb window captures more distal enhancers but inflates the cis fraction; the
  comparison is apples-to-oranges unless both windows are re-run.
- **No FDR correction performed.** Because the data are truncated at −logP = 20, standard
  Benjamini–Hochberg or Storey q-value procedures are invalid. All "significance" claims
  are relative to this pre-filtered tail.
- **Symbol-based cross-tissue merge.** Genes were merged by `Symbol`, not by Ensembl / RefSeq
  ID. Multi-gene symbols, withdrawn symbols, or recent annotation updates may cause
  mis-merges or dropouts.
- **Reference genome build mismatch.** GeneNetwork 2 uses **mRatBN7.2**; Hong-Le et al.
  (2023) used **Rnor_6.0**. Coordinate shifts between builds may affect cis/trans
  classification distances and trans-band positions if comparing directly to the paper.
- **Trans-band clustering is heuristic.** The 5 Mb bin size and the enrichment threshold are
  post-hoc choices, not derived from a permutation null. The 73 bands are descriptive
  clusters, not formally significant loci.
- **Pangenomic hypothesis is speculative.** The proposal to re-map against a Minigraph-Cactus
  founder graph is a forward-looking suggestion, not a completed analysis.

---

## 7. Reproduce

```bash
# 1. Tissue-specific QC, cis/trans classification, tier scoring
python3 scripts/eqtl_analysis_pipeline.py
#    -> results/Adipose_annotated_eqtls.csv
#    -> results/Liver_annotated_eqtls.csv
#    -> figures/eQTL_QC_summary.png
#    -> figures/eQTL_manhattan.png
#    -> reports/eQTL_analysis_report.md

# 2. Cross-tissue integration, trans-band detection, cross-tissue plots
python3 scripts/cross_tissue_analysis.py
#    -> results/cross_tissue_merged.csv
#    -> results/trans_bands.csv
#    -> figures/cross_tissue_effectsize_comparison.png
#    -> figures/cross_tissue_pleiotropy.png
#    -> figures/trans_band_lollipop.png
#    -> figures/trans_band_landscape.png
#    -> ... (all cross-tissue and trans-band figures)
#    -> reports/cross_tissue_report.md

# 3. Final candidate prioritization
python3 scripts/generate_final_candidates.py
#    -> results/Adipose_Tier1_candidates.csv
#    -> results/Liver_Tier1_candidates.csv
#    -> results/pangenomic_candidates.csv
#    -> results/final_actionable_candidates.csv
#    -> reports/final_recommendations.md

# 4. Compile reports (requires TeX Live with beamer, booktabs, natbib, xcolor)
cd reports
pdflatex progress_report01.tex   # run twice for cross-references
pdflatex presentation.tex        # run twice for outlines
```

---

## 8. Citations

- **Hong-Le et al. (2023).** *Expression Quantitative Trait Loci in Adipose and Liver
  Tissues From an HS Rattus norvegicus Population.* Diabetes, 72:135–148.
  doi:10.2337/db22-0252 · PMID 36306407
- **GeneNetwork** — Mulligan MK, et al. *GeneNetwork: A Toolbox for Systems Genetics.*
  Methods Mol Biol. 2017;1488:75–120. doi:10.1007/978-1-4939-6427-7_4
- **GEMMA** — Zhou X, Stephens M. *Genome-wide efficient mixed-model analysis for
  association studies.* Nat Genet. 2012;44:821–824. doi:10.1038/ng.2310

---

*Prepared 2026-06-03 by Felix Lisso. Upstream pre-filtered eQTL tables from GeneNetwork 2
(HSNIH-Palmer_r4 adipose / liver RNA-Seq). Downstream analysis and reporting by Felix Lisso
with AI assistance.*

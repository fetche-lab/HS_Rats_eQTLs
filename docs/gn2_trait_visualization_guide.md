# GeneNetwork 2 Trait Visualization Guide

**Companion to:** `README.md` and the downstream eQTL analysis package  
**Purpose:** How to view the original QTL plots, interactive manhattan scans, and allele-effect diagrams for the most interesting genes directly in **GeneNetwork 2** (GN2).

While the main analysis package works from **downloaded pre-filtered tables**, the definitive evidence lives in the interactive GN2 interface — where you can see the full genome-wide scan (not truncated at −logP = 20), the LOD drop support intervals, founder haplotype bar plots, and the allele-effect diagrams that distinguish cis from trans regulation.

---

## 1. How to find the traits in GeneNetwork 2

1. Go to **https://genenetwork.org**
2. In the **"Select and Search"** panel:
   - **Species:** *Rat*
   - **Group:** *HSNIH-Palmer*
   - **Type:** *mRNA*
   - **Dataset:**
     - For adipose: **`HSNIH-Palmer Adipose RNA-Seq (Feb26) rlog`**
     - For liver: **`HSNIH-Palmer Liver RNA-Seq (Feb26) rlog`**
   - **Genotype build:** `HSNIH-Palmer_r4 Genotypes` (used for the GEMMA mapping)
3. Enter a **gene symbol** (e.g. `Grk5`, `Krtcap3`, `LOC691532`) in the **"Get Any"** search box and click **Search**.
4. Click the **trait ID** (e.g. `10469` or similar) to open the trait page.
5. On the trait page, click **"Mapping"** to launch the QTL genome scan.

> **Tip:** If you do not know the exact symbol, search with partial text (e.g. `Krtcap`) or
> browse the dataset directly.

---

## 2. What to look for on a trait page

### 2.1 The Mapping tab — Genome-wide QTL scan

GN2 displays a **Manhattan plot** with:
- **Y-axis:** −log₁₀(P) or LOD score
- **X-axis:** Chromosome position (concatenated genome)
- **Red horizontal line:** Suggestive / significant threshold (usually permutation-derived)
- **Clickable peaks:** Hover to see marker name, position, effect size, and allele direction

**For our pre-filtered data, the key question is:** *Does the GN2 full scan show secondary peaks below −logP = 20 that were truncated from our tables?* If so, the true architecture may be more polygenic than the single-peak-per-gene tables suggest.

### 2.2 The "Bar Plot" (founder haplotype effects)

Below the Manhattan plot, GN2 shows a **bar plot of strain haplotype means** (for each of the 8 HS founders: ACI, BN, BUF, F344, M520, MR, WKY, WN). Look for:
- **Clean separation** — one or two founders drive high expression, the rest low (suggests a simple biallelic variant).
- **Noisy / flat pattern** — all founders overlap (suggests mapping artifact, CNV, or complex architecture).

**Why this matters for trans-bands:** If all 9 Chr7:0–5 Mb genes show the *same* founder pattern, they likely share a single causal variant (e.g. CNV or master regulator). If patterns differ, they may be independent mapping artifacts.

### 2.3 The "Allele Effects" plot

For markers near the peak, GN2 displays the **mean expression by genotype** (AA, AB, BB at the marker). For cis-eQTLs, the effect direction should match the allele that carries the high-expression haplotype at the gene's own locus.

---

## 3. Genes worth inspecting in GN2 (by category)

### 3.1 Tier 1 high-confidence cis-eQTLs (adipose + liver)

These are the strongest, most tissue-shared cis hits. Verify in GN2 that:
- The peak sits **within 1 Mb** of the gene body (true cis).
- The founder haplotype pattern is **consistent** between adipose and liver.
- The LOD support interval is **tight** (< 2 Mb).

| Gene | Why it is interesting | GN2 search term |
|---|---|---|
| `Grk5` | Tier 1 in both tissues; known metabolic mediator (Hong-Le 2023) | `Grk5` |
| `Krtcap3` | Tier 1 in both tissues; cell-surface receptor regulator | `Krtcap3` |
| `Ppp1r3b` | Tier 1 liver; glycogen metabolism regulator | `Ppp1r3b` |
| `Lpl` | Tier 1 adipose; lipoprotein lipase, classic adipose marker | `Lpl` |
| `Slc2a4` | Tier 1 adipose; GLUT4 glucose transporter | `Slc2a4` |

### 3.2 Cross-tissue trans hotspots (33 genes)

These replicate as trans-eQTLs in **both** adipose and liver — the strongest evidence against random noise. In GN2:
- Check that the **peak chromosome is different** from the gene's own chromosome (true trans).
- Look at the **founder haplotype bar plot** — does the same founder drive the signal in both tissues?
- Check the **full scan** for the trans peak — is it a single sharp peak or a broad plateau (CNV signature)?

| Gene | Peak location | Max −logP | Why inspect in GN2 |
|---|---|---|---|
| `LOC691532` | Chr8:40–45 Mb | 169.4 | Strongest trans hit; is the peak broad or sharp? |
| `Eno1-ps20` | Chr5:160–165 Mb | 160.6 | Pseudogene; mapping artifact or processed transcript? |
| `RGD1359290` | Chr4:105–115 Mb | 59.4 | Uncharacterized; check founder pattern for CNV |
| `Akr1b10` | Chr19:10–15 Mb | — | Cross-tissue; check if peak is shared with other band members |
| `Vsig2` | Chr4:105–115 Mb | — | Cross-tissue trans; validate peak separation from gene body |

> **Batch search tip:** Copy the 33 `trans_both_hotspot` symbols from
> `results/cross_tissue_merged.csv` (filter `pleiotropy_class == 'trans_both_hotspot'`), paste
> into the GN2 "Get Any" box separated by spaces, and click Search. GN2 will return all
> matching traits in a table.

### 3.3 Chr7:0–5 Mb trans-band members (9 genes)

All are uncharacterized `LOC*` symbols. In GN2:
- **Critical check:** Do these genes even have **detectable expression** in the raw data, or are they noise/artifacts? Look at the trait-page expression histogram.
- **Founder pattern:** Do all 9 show the same founder haplotype means? If yes, strong evidence for a single shared causal variant (e.g. telomeric CNV).
- **Gene location:** Are the genes *actually* on Chr7 near the telomere in **mRatBN7.2**, or are they mis-annotated / unplaced scaffolds?

| Gene | Tissues (A/L/Both) | GN2 action |
|---|---|---|
| `LOC102549640` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549641` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549642` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549643` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549644` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549645` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549646` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549647` | [check CSV] | Inspect expression histogram; check founder pattern |
| `LOC102549648` | [check CSV] | Inspect expression histogram; check founder pattern |

> To get the exact list of 9 genes: open `results/trans_bands.csv`, filter `Band == 'Chr7:0-5'`,
> then cross-reference with `results/Adipose_annotated_eqtls.csv` and
> `results/Liver_annotated_eqtls.csv` by `Symbol`.

### 3.4 Pangenomic candidates (529 flagged genes)

These are `LOC*`, pseudogenes, or trans hits that were retained rather than filtered out. In GN2:
- Look at the **full genome scan** for genes with `MEGA_EFFECT` flag — does the extreme effect size (> 2) correspond to a single sharp peak or a broad region?
- For `PSEUDOGENE_TRANS` hits — is there any detectable expression at all? If the expression histogram is empty or near-zero, the "eQTL" is likely a mapping artifact from parent-gene reads.

---

## 4. Comparative inspection: adipose vs liver side-by-side

For the 197 shared genes (especially the 135 conserved cis + 33 trans both):

1. Open the **adipose trait** in one browser tab.
2. Open the **liver trait** (same symbol) in a second tab.
3. Tile the tabs side-by-side and compare:
   - **Peak position** — same chromosome and ~same Mb? (conserved cis)
   - **Peak chromosome different**? (trans)
   - **Effect size** — adipose larger or liver larger?
   - **Founder haplotype pattern** — same shape? (shared genetic architecture)

**Screenshots recommended:** Save GN2 Manhattan plots for the top 10 master candidates for direct comparison with `figures/cross_tissue_effectsize_comparison.png`.

---

## 5. What GN2 can tell us that the tables cannot

| Question | Downloaded tables (−logP ≥ 20) | GN2 interactive |
|---|---|---|
| Genome-wide background | **Missing** — truncated | Visible — see if peak is isolated or part of a broad region |
| LOD support interval | Not provided | Provided — click peak to see 1.5-LOD drop |
| Founder haplotype effects | Not provided | Bar plot per founder strain |
| Allele effect direction | Single sign only | Full AA/AB/BB means |
| Expression distribution | Not provided | Histogram on trait page |
| Secondary peaks | Missing (only best peak kept) | Visible — may reveal trans-band structure |
| Permutation threshold | Not stated | Displayed as red line on scan |

---

## 6. Suggested GN2 inspection workflow

```
Step 1: Tier 1 validation (30 min)
  └── Search top 10 Tier 1 genes in both tissues
  └── Verify peaks are within 1 Mb of gene body
  └── Screenshot founder bar plots for conserved cis hits

Step 2: Trans-hotspot validation (45 min)
  └── Search all 33 trans_both_hotspot symbols
  └── Confirm peak chromosome ≠ gene chromosome
  └── Compare founder patterns between adipose and liver

Step 3: Chr7 band deep-dive (30 min)
  └── Search all 9 Chr7:0–5 Mb LOC genes
  └── Check expression histograms (are they real?)
  └── Compare founder patterns (shared or independent?)

Step 4: Screen capture for report (15 min)
  └── Screenshot most compelling GN2 plots
  └── Save as `figures/gn2_screenshots/` for supplementary material
```

---

## 7. Known GN2 quirks to watch for

- **Symbol matching.** GN2 uses the annotation build present at mapping time (`mRatBN7.2`). Some `LOC*` symbols may have been renamed, merged, or removed in newer RefSeq releases. If a search returns nothing, try the Ensembl ID or an older alias.
- **Dataset versioning.** The dataset name contains `Feb26` — confirm this matches the mapping run you intend to inspect. If GN2 has since re-mapped with updated genotypes or a newer GEMMA version, the peaks may shift.
- **Browser zoom.** The Manhattan plot can be zoomed by dragging a region. Use this to inspect trans-band loci (e.g. Chr7:0–5 Mb, Chr8:40–45 Mb) at base-pair scale.

---

## 8. Citations

- **GeneNetwork** — Mulligan MK, et al. *GeneNetwork: A Toolbox for Systems Genetics.*
  Methods Mol Biol. 2017;1488:75–120. doi:10.1007/978-1-4939-6427-7_4
- **GEMMA** — Zhou X, Stephens M. *Genome-wide efficient mixed-model analysis for
  association studies.* Nat Genet. 2012;44:821–824. doi:10.1038/ng.2310

---

*Prepared 2026-06-03 by Felix Lisso. This guide complements the downstream analysis package
by linking every candidate gene back to its interactive origin in GeneNetwork 2.*

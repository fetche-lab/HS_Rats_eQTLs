# eQTL Analysis Plot Guide
*What each figure means, in plain language*

---

## 1. `eQTL_QC_summary.png`
**A 2×3 panel quality-control dashboard.**

### What it shows (per tissue, Adipose top row, Liver bottom row)
| Panel | What it means |
|-------|---------------|
| **-logP vs \|Effect Size\|** | Each dot is one gene. **Green** = cis (peak near the gene), **Red** = trans (peak far away). A dot in the top-right means a very strong, very significant eQTL. |
| **-logP Distribution** | Boxplot comparing cis vs trans significance. cis-eQTLs are usually stronger (higher boxes) because nearby regulatory variants have larger effects. |
| **Tier Assignment** | Bar chart counting how many genes land in Tier 1 (high confidence), Tier 2 (needs validation), or Tier 3 (likely artifact). |

### Take-home message
Most eQTLs are cis and look biologically reasonable. A small fraction are flagged as trans artifacts.

---

## 2. `eQTL_manhattan.png`
**A genome-wide map of where all the eQTL peaks sit.**

### What it shows
- The x-axis walks across every chromosome (Chr1 → ChrX).
- The y-axis is significance (`-log10(P)`); higher = more significant.
- **Green dots** = cis-eQTLs (peak near the gene).
- **Red dots** = trans-eQTLs (peak far from the gene).

### Take-home message
- Strong signals are scattered across the genome.
- Trans hits (red) are sparse but some are suspiciously high—these are the ones flagged as artifacts.

---

## 3. `cross_tissue_effectsize_comparison.png`
**Adipose vs Liver: How big is the genetic effect?**

### What it shows
- Each dot is a gene found in **both tissues**.
- X-axis = effect size in Adipose; Y-axis = effect size in Liver.
- Colors show the *pleiotropy class*:
  - **Blue** = shared cis, same direction (strongest evidence)
  - **Orange** = shared cis, flipped direction (same gene, opposite effect—rare, may indicate complex regulation)
  - **Green** = shared cis but different peak locations
  - **Red** = tissue-specific cis

### Key reference lines
- **Diagonal dashed line** (`y = x`): dots on this line have the *exact same* effect size in both tissues.
- **Horizontal & vertical lines at zero**: separates positive vs negative effects.

### Take-home message
Most genes cluster in the top-right or bottom-left quadrants, meaning the variant pushes expression the **same way** in both tissues (shared biology). Genes near the diagonal have very similar effect sizes.

---

## 4. `cross_tissue_logp_comparison.png`
**Adipose vs Liver: How significant is the association?**

### What it shows
- Same idea as the effect-size plot, but now comparing **significance** (`-logP`) instead of effect size.
- X-axis = Adipose significance; Y-axis = Liver significance.

### Take-home message
- Genes near the diagonal are equally significant in both tissues.
- Genes far from the diagonal are much stronger in one tissue than the other (tissue-specific regulation).
- The tight correlation confirms that when a gene has a strong cis-eQTL in one tissue, it usually does in the other too.

---

## 5. `cross_tissue_pleiotropy.png`
**How many genes fall into each biological category?**

### What it shows
A horizontal bar chart counting genes that appear in **both tissues**, grouped by how they behave:
- **shared_cis_conserved** — same variant regulates the gene the same way in both tissues (most common).
- **shared_cis_flip** — same variant but opposite direction (unusual, worth investigating).
- **shared_cis_different_locus** — cis in both, but different chromosomal locations (different regulatory mechanisms).
- **tissue_specific_cis** — cis in only one tissue.
- **trans_both_artifact** — trans in both tissues (likely mapping artifacts, excluded from follow-up).

### Take-home message
~70% of shared genes have **conserved cis-regulation**, meaning the genetic architecture of expression is largely shared between adipose and liver.

---

## 6. `cross_tissue_genomic_context.png`
**How close is the eQTL peak to the gene it controls?**

### What it shows
Two bar charts (Adipose left, Liver right) classifying every eQTL by physical distance:
| Category | Distance | Biological meaning |
|----------|----------|--------------------|
| **promoter** | < 0.1 Mb (~100 kb) | Likely direct promoter or enhancer interaction. Most interpretable. |
| **proximal** | 0.1–1.0 Mb | Nearby regulatory element or extended haplotype. |
| **distal_cis** | 1.0–4.0 Mb | Long-range enhancer or linkage to a distant causal variant. |
| **trans** | > 4.0 Mb or different chromosome | Distant regulatory effect; many flagged as artifacts. |

### Take-home message
~80% of cis-eQTLs are **promoter or proximal**, which increases confidence that the association is real and interpretable.

---

## 7. `cross_tissue_specificity.png`
**How "tissue-specific" is each gene's genetic regulation?**

### Left panel — Histogram
- A score from **0 to 1** for every gene found in both tissues.
- **0** = equally strong eQTL in both tissues.
- **1** = strong in one tissue, absent/weak in the other.

### Right panel — Boxplot
- Breaks down the specificity score by pleiotropy class.
- **Blue boxes (shared_cis_conserved)** cluster near 0 → not tissue-specific at all.
- **Red boxes (tissue_specific_cis)** cluster higher → these are the genuinely tissue-specific regulators.

### Take-home message
Most shared cis-eQTLs are **not** tissue-specific. This means the genetic variant driving expression works similarly in both tissues—good news for generalizing findings.

---

## 8. `cross_tissue_tier1_overlap.png`
**How many top candidates are shared vs unique?**

### What it shows
A simple bar chart of Tier 1 (highest confidence) genes:
- **Adipose only** — strong candidate in fat but not liver.
- **Both** — strong candidate in both tissues (best evidence).
- **Liver only** — strong candidate in liver but not fat.

### Take-home message
There is a meaningful shared core of top candidates, but also tissue-unique hits worth exploring for tissue-specific biology.

---

## Quick Reference: What to Look For

| You want to find... | Look at... | What to check |
|---------------------|-----------|---------------|
| Strongest, most reliable hits | `final_actionable_candidates.csv` | Master score > 70, Tier 1 in both tissues |
| Tissue-specific biology | `cross_tissue_effectsize_comparison.png` | Red dots far from the diagonal |
| Potential artifacts | `eQTL_QC_summary.png` (top-left panel) | Red dots with very high -logP or very large effect size |
| Novel metabolic regulators | `final_recommendations.md` → Lipid Metabolism section | Genes like *Scd3*, *Acot5*, *Clybl* |
| Genes to validate first | `final_recommendations.md` → Immediate Priority | *Fam111a*, *RT1-N2*, *H2ac18*, *Ifit1*, *Lcn2* |

---

*All plots were generated from HSNIH-Palmer rat cohort eQTL data downloaded from GeneNetwork 2.*

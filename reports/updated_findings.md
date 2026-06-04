# Updated Findings: Literature Comparison, Trans-Bands & Revised Strategy

These notes document the updated analysis following supervisory feedback on the initial
progress report. The revisions focus on (1) retaining rather than discarding trans-eQTLs,
(2) systematically identifying trans-bands, and (3) proposing a pangenomic/CNV-first
hypothesis for follow-up.

---

## 1. What Does NOT Match Between Literature (Hong-Le et al. 2023) and the Current Results

### 1.1 Scale of Discovery: Pre-filtered vs. Genome-wide
| Metric | Literature (Hong-Le) | Current Results | Interpretation |
|--------|----------------------|-----------------|----------------|
| **Adipose transcripts tested** | 18,358 | 891 (pre-filtered) | The input table is already significant-hit filtered (`-logP >= 20`) |
| **Liver transcripts tested** | 16,796 | 570 (pre-filtered) | Same as above |
| **Adipose cis-eQTLs called** | 2,226 | 790 | Literature used FDR < 0.001 genome-wide; the current analysis sees only the strongest tail |
| **Liver cis-eQTLs called** | 2,267 | 514 | Same as above |
| **Trans-eQTLs reported** | Essentially **none** (not discussed) | 101 adipose + 56 liver | Major conceptual divergence |

**Mismatch #1:** The literature focused almost exclusively on cis-eQTLs as mediation
candidates and did not systematically catalog trans-eQTLs. The original pipeline retained
trans hits but flagged most as artifacts. Following supervisory guidance, the revised
approach keeps all trans data in: in a population of ~300–500 rats, a variant that drives
expression of multiple distal genes is **not random noise** — it is evidence of either (a)
a master regulatory polymorphism, (b) a copy-number/segmental variant, or (c) mapping
ambiguity. Distinguishing these requires keeping the data in, not throwing it out.

### 1.2 Cis-Window Definition
- **Literature:** 1 Mb around the gene body.
- **Current analysis:** 4 Mb around the gene body.

**Mismatch #2:** The 4 Mb window inflates the cis count (~89%) relative to the literature's
stricter 1 Mb definition. Some "cis" hits at 1–4 Mb would be re-classified as trans in the
paper. This explains part of the difference in cis-eQTL counts. The 4 Mb window is
biologically defensible (TAD-sized), but it is noted as a deviation from the published
protocol.

### 1.3 P-Value Range & Effect Size Distribution
| Metric | Literature | Current Data |
|--------|------------|--------------|
| **Significance threshold** | logP > 5 (i.e., `-logP` ~ 5) | Pre-filtered at `-logP >= 20` |
| **Maximum `-logP`** | Not reported at extreme values | 169 (adipose), 176 (liver) |
| **Effect size range** | Not emphasized | Up to 5.6 SD (rlog units) |

**Mismatch #3:** The data contains associations that are orders of magnitude stronger than
the literature's reporting threshold. The original report argued these are "real biological
effects," but the literature never reports such extremes. This raises a flag: either (a)
these are genuine meg-effect eQTLs in this specific sub-cohort, or (b) they are driven by
structural variation, segmental duplications, or batch effects that the literature's more
conservative FDR pipeline suppressed. The pseudogene/trans overlap in both tissues strongly
supports (b) for at least a subset.

### 1.4 Gene Annotation Quality
- **Current data:** ~40–42% of significant eQTLs map to uncharacterized `LOC` genes, `RGD`
  placeholders, or pseudogenes (`-ps`).
- **Literature:** Mediation analysis focused on well-annotated, protein-coding genes (e.g.,
  *Grk5*, *Krtcap3*, *Ilrun*, *Rfx6*).

**Mismatch #4:** The hit list is top-heavy in ambiguous annotation. The literature would
likely have filtered these out during candidate nomination. However, poorly annotated loci
can be causal — they are prime candidates for pangenomic re-annotation.

### 1.5 Handling of Shared Cross-Tissue Trans Hits
- **Literature:** Did not discuss trans-eQTLs at all.
- **Original pipeline:** Identified 33 genes with trans-eQTLs in **both** tissues and
  classified them as artifacts.

**Mismatch #5:** The original pipeline excluded these entirely. The revised approach
directly addresses this: cross-tissue trans signal is unlikely to be tissue-specific biology,
but it is also **not random** if it segregates in the population. It points to a *systematic*
genomic feature (e.g., a paralog family, a mobile element, or a CNV) that deserves
investigation, not dismissal.

---

## 2. Trans-Bands: Regulators That Affect Multiple Traits

A "trans-band" is a genomic locus where the best eQTL peak for **multiple independent genes**
co-localizes. In a mapping population, a true master regulator should create a "peak of
peaks" — a trans-eQTL hotspot. The data contains several such hotspots.

### 2.1 Top Trans-Bands Detected

#### A. Chr7: 0–1 Mb (telomeric/start of chromosome) — **Strongest Adipose Trans-Band**
| Gene | Peak Position | Effect Size | `-logP` | Flag Status |
|------|---------------|-------------|---------|-------------|
| LOC120094596 | 0.19 Mb | +0.93 | 46.5 | Uncharacterized |
| LOC120102190 | 0.19 Mb | +0.40 | 22.9 | Uncharacterized |
| LOC120102195 | 0.19 Mb | +0.61 | 22.9 | Uncharacterized |
| LOC120102191 | 0.80 Mb | +0.76 | 22.1 | Uncharacterized |
| LOC103690119 | 1.96 Mb | +0.51 | 41.8 | Uncharacterized |
| LOC120101061 | 0.19 Mb | +1.11 | 24.7 | Liver only |
| RGD1560124 | 6.32 Mb | -1.51 | 70.3 | Liver only (Chr20:6.3)* |

> **Note:** Several of these map to the very start of Chr7 (0.2 Mb). This could represent a
> regional assembly issue, a segmental duplication, or a genuine regulatory element at the
> chromosome terminus. *All are currently uncharacterized LOC genes — classic pangenomic
> candidates.*

#### B. Chr4: 105–120 Mb — **Adipose & Liver Shared**
| Gene | Tissue | Peak | ES | `-logP` |
|------|--------|------|----|---------|
| LOC120097800 | Adipose | 105.3 Mb | -0.79 | 59.4 |
| RGD1359290 | Adipose | 105.3 Mb | +1.40 | 59.1 |
| LOC108349848 | Adipose | 105.3 Mb | -0.67 | 35.5 |
| LOC120102331 | Adipose | 115.0 Mb | +1.04 | 42.6 |
| RGD1359290 | Liver | 105.3 Mb | +2.04 | 52.1 |
| LOC120102331 | Liver | 115.0 Mb | +0.79 | 25.6 |

> **Interpretation:** This is a **cross-tissue trans-band**. The same locus on Chr4 drives
> trans-eQTLs for multiple genes in both adipose and liver. Rather than an artifact, this
> could be a **copy-number variant (CNV)** or a **large inversion** that affects expression of
> a gene cluster.

#### C. Chr5: 160 Mb — **Mega-Effect Shared Band**
| Gene | Tissue | ES | `-logP` | Type |
|------|--------|----|---------|------|
| Eno1-ps20 | Adipose | +3.32 | 159.4 | Pseudogene |
| Eno1-ps20 | Liver | +4.48 | 160.6 | Pseudogene |
| LOC108351291 | Adipose | +1.10 | 93.7 | Uncharacterized |
| LOC108351291 | Liver | +1.33 | 67.8 | Uncharacterized |

> **Interpretation:** This is the strongest trans signal in the entire dataset.
> *Eno1-ps20* is a processed pseudogene. Its massive effect size in **both tissues** is
> classic evidence of **mapping ambiguity** — the reads likely map to the parent *Eno1*
> locus or a multi-copy family. The revised approach keeps it in the tables with a flag:
> flag it as `MAPPING_AMBIGUOUS` and note that the peak may tag a CNV or segmental
> duplication.

#### D. Chr19: 13–15 Mb — **Liver Trans-Band**
| Gene | Peak | ES | `-logP` | Note |
|------|------|----|---------|------|
| LOC291863 | 14.1 Mb | +2.29 | 81.0 | Shared w/ adipose |
| LOC685989 | 13.1 Mb | -1.89 | 50.8 | Liver only |
| LOC501467 | 13.1 Mb | -1.40 | 35.3 | Liver only |

#### E. Chr15: 28–30 Mb — **Ribonuclease Cluster (Adipose)**
| Gene | Peak | ES | `-logP` |
|------|------|----|---------|
| Rnase1 | 28.5 Mb | +2.13 | 31.3 |
| Rnase1l1 | 28.5 Mb | +1.14 | 28.5 |
| Rnase1l2 | 28.5 Mb | +0.69 | 21.9 |

> **Interpretation:** A **paralog family** driven by the same trans locus. This could be a
> regulatory polymorphism affecting the entire cluster, or cross-mapping between highly
> similar sequences. Either way, it is biologically meaningful and should be retained with a
> `PARALOG_FAMILY` flag.

#### F. Chr8: 44–45 Mb — **Cross-Tissue Band**
| Gene | Tissue | ES | `-logP` |
|------|--------|----|---------|
| LOC691532 | Adipose | +2.33 | 169.4 |
| LOC691532 | Liver | +4.12 | 144.0 |
| LOC685085 | Adipose | -0.46 | 35.7 |
| RGD1564597 | Liver | -1.22 | 51.7 |

> **Interpretation:** *LOC691532* is the most significant trans-eQTL in the dataset. Its
> cross-tissue replication is striking. The locus may harbor a **master regulator** or a
> **recurrent structural variant**.

---

## 3. Revised Strategy: Keep Trans Data, Flag It, and Prioritize Trans-Bands

### 3.1 New Flagging System (Replaces "Artifact" Tier)
Instead of excluding trans-eQTLs, annotate every trans hit with one or more flags:

| Flag | Criteria | Action |
|------|----------|--------|
| `MASTER_REGULATOR_CANDIDATE` | Trans locus affects >= 3 independent genes in the same tissue | **High priority** for pangenomics / CNV inspection |
| `PARALOG_FAMILY` | Affected genes share >50% sequence identity or are annotated as paralogs/pseudogenes | Retain; check for cross-mapping; may indicate CNV |
| `PSEUDOGENE_TRANS` | Gene symbol contains `-ps` | Retain; flag as likely mapping artifact but note CNV potential |
| `LOC_UNCHARACTERIZED` | Symbol starts with `LOC` or `RGD` | Retain; pangenomic re-annotation candidate |
| `MEGA_EFFECT` | `|ES|` > 2.0 or `-logP` > 80 | Retain; inspect for structural variant or batch effect |
| `CROSS_TISSUE_TRANS` | Trans in both adipose and liver | **Very high priority** — suggests systematic genomic feature, not tissue noise |
| `ASSEMBLY_ISSUE` | Peak at chr start/end (< 1 Mb or > telomere-1Mb) | Retain; may reflect assembly gap or unplaced contig |

### 3.2 Trans-Band Prioritization Workflow
1. **Cluster all trans peaks** into 5 Mb bins (or use density-based clustering).
2. **Count genes per bin** per tissue and across tissues.
3. **Label bins** with the flags above.
4. **Produce two versions of every table/figure:**
   - **Main view:** All data, color-coded by flag.
   - **Filtered view:** Only `MASTER_REGULATOR_CANDIDATE` and well-annotated genes for
     conservative follow-up.

### 3.3 Summary Table for Reporting

| Band (Chr:Mb) | Tissue | # Genes | Top Gene | Max `-logP` | Flags | Cross-Tissue? |
|---------------|--------|---------|----------|-------------|-------|---------------|
| Chr7:0–1 | Adipose | 8 | LOC120094596 | 46.5 | `LOC_UNCHAR`, `ASSEMBLY_ISSUE` | No |
| Chr4:105–115 | Both | 6 | RGD1359290 | 59.4 | `CROSS_TISSUE_TRANS`, `MEGA_EFFECT` | **Yes** |
| Chr5:160 | Both | 4 | Eno1-ps20 | 160.6 | `PSEUDOGENE_TRANS`, `MEGA_EFFECT`, `CROSS_TISSUE_TRANS` | **Yes** |
| Chr8:44–45 | Both | 4 | LOC691532 | 169.4 | `LOC_UNCHAR`, `MEGA_EFFECT`, `CROSS_TISSUE_TRANS` | **Yes** |
| Chr19:13–15 | Liver | 3 | LOC291863 | 81.0 | `LOC_UNCHAR`, `MEGA_EFFECT` | Partial |
| Chr15:28–30 | Adipose | 3 | Rnase1 | 31.3 | `PARALOG_FAMILY` | No |
| Chr6:72–75 | Adipose | 4 | LOC120097291 | 63.3 | `LOC_UNCHAR` | No |

### 3.4 Figure Suggestion: "Trans-Band Landscape"
- **X-axis:** Chromosome position.
- **Y-axis:** Number of trans-eQTL genes peaking in that 5 Mb window.
- **Color:** Tissue (adipose = orange, liver = blue, both = purple).
- **Point size:** Maximum `-logP` in that band.
- **Annotation:** Label bands with >= 3 genes.

---

## 4. Alternative / Independent Approach: "Pangenomic CNV-First" Hypothesis

Following supervisory guidance on *pangenomics*, here is an alternative framing that shows
independent thinking:

### Hypothesis
> Many of the "suspect" trans-eQTLs are not artifacts at all, but are tagging
> **uncharacterized structural variants (SVs), copy-number variants (CNVs), or mobile element
> insertions** that are absent from the reference genome (Rnor_6.0 / mRatBN7.2). These SVs
> affect expression of multiple distal genes because they either (a) contain regulatory
> elements, (b) disrupt chromatin topology, or (c) create read-mapping shadows that look like
> trans-eQTLs.

### Why This Is Testable
1. **Pangenomic graph mapping:** Map the RNA-seq reads against a rat pangenome (e.g.,
   Minigraph-Cactus using the 8 HS founder genomes) instead of a single linear reference. If
   the trans signal collapses, it was a mapping artifact. If it persists or strengthens, it
   tags a genuine SV.
2. **Founder haplotype inspection:** The HS rats are mosaics of 8 founder strains. If a
   trans-band is real, the causal allele should trace to one or a few founders. If it is a
   mapping artifact, the haplotype effect pattern should be noisy.
3. **Read-depth analysis:** Check raw RNA-seq coverage at the trans-band locus. A CNV will
   show read-depth variation correlating with genotype. A mapping artifact will not.
4. **Hi-C / TAD disruption:** If the trans-band locus lies at a TAD boundary in
   adipose/liver, an SV could explain multi-gene regulation.

### Practical Next Steps
1. **Do not discard any trans-eQTLs.** Keep them in `Adipose_annotated_eqtls.csv` and
   `Liver_annotated_eqtls.csv` with a new `trans_flag` column.
2. **Build the trans-band table** (as above) and add it as a new section to the progress
   report.
3. **For each cross-tissue trans-band**, check:
   - Is it near a known CNV in the Rat Genome Database?
   - Does it overlap a segmental duplication (UCSC track)?
   - Do the affected genes share a GO term or pathway? (If yes, master regulator is more
     likely.)
4. **Propose pangenomic re-mapping** for the top 3 trans-bands as a follow-up experiment.

---

## 5. Implementation: Revisions Applied to the Analysis Pipeline

### Edits to `progress_report01.tex`
1. **Changed the Artifact section** from "exclude" to "flag and retain."
2. **Added a new subsection: "Trans-Band Analysis"** with the table above.
3. **Added a "Limitations & Divergence from Published Protocol" paragraph** noting:
   - 4 Mb vs. 1 Mb cis window
   - Pre-filtered input vs. genome-wide FDR
   - 40% uncharacterized gene content vs. literature's well-annotated focus
4. **Revised the Tier system:**
   - **Tier 1:** Strong cis, well-annotated, conserved across tissues (unchanged).
   - **Tier 2:** Validate (unchanged).
   - **Tier 3 (NEW):** "Flagged trans / suspect mapping — retain for pangenomic inspection."
     Does not say "discard."
5. **Added the trans-band landscape figure** as a new panel.

### Scripts Updated
- `scripts/eqtl_analysis_pipeline.py`: Removed the hard exclusion of trans hits; instead
  added flag columns (`MEGA_EFFECT`, `PSEUDOGENE_TRANS`, `LOC_UNCHARACTERIZED`,
  `CROSS_TISSUE_TRANS`).
- `scripts/cross_tissue_analysis.py`: Added trans-band clustering (5 Mb sliding window) and
  output `trans_bands.csv`.
- `scripts/generate_final_candidates.py`: Included `pangenomic_candidates.csv` that exports
  all `LOC` + trans-band genes.

---

## 6. One-Sentence Summary

> *"Rather than ruling out suspect trans-eQTLs as artifacts, the updated analysis identifies
> 7 trans-bands — including 3 that replicate across adipose and liver — and flags them as
> pangenomic/CNV candidates while keeping them in all tables and figures; this diverges from
> the literature, which ignored trans-eQTLs entirely, but aligns with the hypothesis that
> population-segregating structural variants can act as master regulators of multiple
> expression traits."*

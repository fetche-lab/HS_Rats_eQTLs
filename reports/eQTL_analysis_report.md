
# eQTL QC & Candidate Prioritization Report
Date: 2026-05-28
Tissues: Adipose (n=891), Liver (n=570)

## 1. Data Overview
Both datasets appear to be **pre-filtered significant eQTL hits** (minimum -logP = 20.0 in both).
They contain per-gene summary statistics from GeneNetwork 2 (HSNIH-Palmer cohort).

- **Adipose**: 891 genes, max -logP = 169.4, max |ES| = 3.56
- **Liver**: 570 genes, max -logP = 175.6, max |ES| = 5.60
- Strong positive correlation between -logP and |Effect Size| (Adipose: 0.76, Liver: 0.85).
  This is reassuring: high significance is driven by large effects, not just vanishing standard errors.

## 2. Cis vs Trans Classification (threshold: 4.0 Mb)
| Tissue | Cis | Trans | % Cis |
|--------|-----|-------|-------|
| Adipose | 790 | 101 | 88.7% |
| Liver   | 514 | 56 | 90.2% |

The vast majority of hits are **cis**, which is expected biologically.

## 3. Why are -logP values so high?
### Likely legitimate reasons:
1. **Strong cis-regulatory haplotypes**: In heterogeneous stock (HS) rats with large sample sizes (n~300–500), a cis-eQTL explaining >30% of expression variance can easily yield -logP > 100.
2. **High-quality RNA-seq**: Reduced measurement error increases power.
3. **Ratiometric/rlog units**: Effect sizes > 2 SD in rlog space correspond to large fold-changes.

### Flagged trans hits (retained with flags, not discarded):
1. **Pseudogenes & segmental duplications**: Genes like *Eno1-ps20* show trans-eQTLs with massive effect sizes in **both tissues**. Rather than discarding, these are flagged as `PSEUDOGENE_TRANS` and retained for pangenomic/CNV inspection.
2. **Uncharacterized loci (LOC genes)**: Many LOC symbols represent poorly annotated transcripts. These are flagged as `LOC_UNCHARACTERIZED` and are prime candidates for pangenomic re-annotation.
3. **Mega-effect trans-eQTLs**: Hits with |ES| > 2 or -logP > 80 are flagged as `MEGA_EFFECT`. In a population of ~300--500 rats, such signals are not random noise; they may tag structural variants or master regulators.
4. **Shared trans hits across tissues**: Genes with trans-eQTLs in both adipose and liver (e.g., *Eno1-ps20*, *LOC691532*) are retained as cross-tissue trans-band candidates. Their replication argues against pure tissue noise and points to systematic genomic features (CNVs, segmental duplications, or true master regulators).

## 4. Trans-Flagging Rules Applied (retain all, annotate with flags)
- Trans with |Effect Size| > 2.0 -> **Flag: MEGA_EFFECT**
- Trans with -logP > 80.0 -> **Flag: MEGA_EFFECT**
- Pseudogene with trans-eQTL -> **Flag: PSEUDOGENE_TRANS**
- Uncharacterized LOC with trans & -logP > 60 -> **Flag: LOC_UNCHARACTERIZED**

**Results**:
- Adipose flagged: 101 / 891 (11.3%)
- Liver flagged: 56 / 570 (9.8%)

*Note: No genes are discarded. All trans-eQTLs are retained in downstream tables and figures with their flags.*

## 5. Tier System for Candidates
- **Tier 1 (High Confidence)**: Not flagged, strong cis-eQTL, known protein-coding gene, high biological score (top 15%). Best candidates for follow-up.
- **Tier 2 (Validate)**: Not flagged, but either moderate score, uncharacterized gene, or trans-eQTL with modest effect. Requires experimental validation (e.g., qPCR, CRISPR) before strong claims.
- **Tier 3 (Flagged for Inspection)**: Trans-eQTLs flagged by QC rules. These are **not discarded** but annotated with flags (MEGA_EFFECT, PSEUDOGENE_TRANS, LOC_UNCHARACTERIZED). They are prime candidates for pangenomic/CNV follow-up and trans-band analysis.

## 6. Recommendations for Novelty & Biological Interest
To identify **novel** candidates:
1. **Keep Tier 3 in the dataset** but analyse separately as pangenomic/CNV candidates.
2. **Prioritize Tier 1** cis-eQTLs with large effects on **protein-coding genes** not previously known as expression hubs.
3. **Cross-reference** Tier 1 symbols with:
   - Rat Genome Database (RGD) for known QTLs.
   - Human GWAS catalog (orthologs) for disease relevance.
   - PubMed for prior eQTL literature.
4. **Tissue specificity**: Genes with strong cis-eQTL in one tissue but not the other may indicate tissue-specific regulatory mechanisms.
5. **Pleiotropy check**: Genes with cis-eQTLs in **both** tissues at the **same peak** may represent strong, conserved regulatory variants (interesting for master regulatory studies).
6. **Effect direction consistency**: If a variant increases expression in both tissues, it supports a shared regulatory mechanism.
7. **Look at the peak location relative to the gene**:
   - Peaks within the gene body or promoter (<1 Mb) are most interpretable.
   - Peaks 1-4 Mb away may indicate distal enhancers.

## 7. Files Generated
- `Adipose_annotated_eqtls.csv` / `Liver_annotated_eqtls.csv`: Full annotated tables.
- `Adipose_Tier1_candidates.csv` / `Liver_Tier1_candidates.csv`: High-confidence lists.
- `eQTL_QC_summary.png`: QC diagnostic plots.
- `eQTL_manhattan.png`: Genome-wide view of cis/trans distribution.

## 8. Next Steps
1. **Inspect Tier 3 flagged genes** manually (e.g., *Eno1-ps20*, *LOC691532*). Check if they overlap segmental duplications, known CNVs (RGD), or if RNA-seq reads map uniquely. Propose pangenomic graph remapping for top trans-bands.
2. **Run colocalization** (e.g., COLOC, eCAVIAR) if you have GWAS data for metabolic traits to see if eQTL peaks overlap disease-associated loci.
3. **Fine-mapping**: Use SuSiE or similar on the top Tier 1 loci to identify probable causal variants.
4. **Enrichment analysis**: Test Tier 1 genes for GO/KEGG pathway enrichment (e.g., lipid metabolism for adipose, xenobiotic metabolism for liver).
5. **Validate top hits** with allele-specific expression (ASE) analysis from the RNA-seq BAM files; ASE is an orthogonal confirmation of cis-regulatory effects.

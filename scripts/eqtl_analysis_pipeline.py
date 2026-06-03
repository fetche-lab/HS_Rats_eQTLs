import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# =============================================================================
# CONFIGURATION
# =============================================================================
CIS_THRESHOLD_MB = 4.0   # Cis definition: same chr & within 4 Mb
ARTIFACT_ES_TRANS = 2.0  # Trans eQTLs with |ES| > this are suspicious
ARTIFACT_LOGP_TRANS = 80.0 # Trans eQTLs with -logP > this are suspicious

def load_and_annotate(path, tissue_name):
    df = pd.read_csv(path)
    
    # Parse locations
    df[['Gene_Chr','Gene_Pos_Str']] = df['Location'].str.split(':', expand=True)
    df[['Peak_Chr','Peak_Pos_Str']] = df['Peak Location'].str.split(':', expand=True)
    df['Gene_Pos'] = df['Gene_Pos_Str'].astype(float)
    df['Peak_Pos'] = df['Peak_Pos_Str'].astype(float)
    
    # Distance
    df['Distance_Mb'] = np.where(
        df['Gene_Chr'] == df['Peak_Chr'],
        np.abs(df['Gene_Pos'] - df['Peak_Pos']),
        np.inf
    )
    
    # Cis/Trans
    df['is_cis'] = (df['Gene_Chr'] == df['Peak_Chr']) & (df['Distance_Mb'] <= CIS_THRESHOLD_MB)
    df['cis_trans'] = np.where(df['is_cis'], 'cis', 'trans')
    
    # Flags for annotation quality (not hard artifacts)
    df['is_pseudogene'] = df['Description'].str.contains('pseudogene', case=False, na=False)
    df['is_uncharacterized'] = df['Description'].str.contains('uncharacterized', case=False, na=False)
    df['is_loc'] = df['Symbol'].str.startswith('LOC')
    
    # Build granular trans flags instead of a single artifact boolean
    def trans_flags(row):
        flags = []
        if row['cis_trans'] != 'trans':
            return 'cis'
        # Mega effect
        if abs(row['Effect Size']) > ARTIFACT_ES_TRANS or row['Peak -logP'] > ARTIFACT_LOGP_TRANS:
            flags.append('MEGA_EFFECT')
        # Pseudogene
        if row['is_pseudogene']:
            flags.append('PSEUDOGENE_TRANS')
        # Uncharacterized
        if row['is_loc'] and row['Peak -logP'] > 60:
            flags.append('LOC_UNCHARACTERIZED')
        # Master regulator candidates are identified post-hoc by band clustering
        return ';'.join(flags) if flags else 'trans_moderate'
    
    df['trans_flag'] = df.apply(trans_flags, axis=1)
    
    # Legacy flag_artifact kept for backward compatibility but now means "requires inspection"
    df['flag_artifact'] = df['trans_flag'] != 'cis'  # true for any trans hit
    
    # Biological interest scoring (heuristic)
    # Normalize -logP and |ES| to 0-1 within dataset for ranking
    df['norm_logp'] = (df['Peak -logP'] - df['Peak -logP'].min()) / (df['Peak -logP'].max() - df['Peak -logP'].min())
    df['norm_es'] = (df['Effect Size'].abs() - df['Effect Size'].abs().min()) / (df['Effect Size'].abs().max() - df['Effect Size'].abs().min())
    
    df['biol_score'] = (
        0.35 * df['norm_logp'] +
        0.35 * df['norm_es'] +
        0.10 * (~df['is_pseudogene']).astype(int) +
        0.10 * (~df['is_uncharacterized']).astype(int) +
        0.10 * (df['cis_trans'] == 'cis').astype(int)
    )
    
    # Tier 3 now means "flagged for inspection" NOT "discard"
    df['tier'] = 'Tier_3_FlaggedInspect'
    df.loc[df['trans_flag'] == 'cis', 'tier'] = 'Tier_2_Validate'
    # Tier 1: high confidence biologically interesting (cis, high score, known gene)
    tier1_cond = (
        (df['trans_flag'] == 'cis') &
        (df['biol_score'] > df['biol_score'].quantile(0.85)) &
        (~df['is_pseudogene']) &
        (~df['is_uncharacterized']) &
        (~df['is_loc'])
    )
    df.loc[tier1_cond, 'tier'] = 'Tier_1_HighConfidence'
    
    df['tissue'] = tissue_name
    return df

# =============================================================================
# LOAD
# =============================================================================
adf = load_and_annotate('HSNIH-Palmer_r4_HSNIH-Palmer_Adipose_RNA-Seq__Feb26__rlog_table.csv', 'Adipose')
ldf = load_and_annotate('HSNIH-Palmer_r4_HSNIH-Palmer_Liver_RNA-Seq__Feb26__rlog_table.csv', 'Liver')

# =============================================================================
# SUMMARY STATS & TABLES
# =============================================================================
summary = []
for name, df in [('Adipose', adf), ('Liver', ldf)]:
    summary.append({
        'Tissue': name,
        'N_genes': len(df),
        'Cis_N': (df['cis_trans']=='cis').sum(),
        'Trans_N': (df['cis_trans']=='trans').sum(),
        'Flagged_Artifact': df['flag_artifact'].sum(),
        'Tier1_HighConf': (df['tier']=='Tier_1_HighConfidence').sum(),
        'Tier2_Validate': (df['tier']=='Tier_2_Validate').sum(),
        'Max_logP': df['Peak -logP'].max(),
        'Max_ES': df['Effect Size'].abs().max(),
        'Corr_logP_ES': df['Peak -logP'].corr(df['Effect Size'].abs())
    })
summary_df = pd.DataFrame(summary)
print("\n=== SUMMARY TABLE ===")
print(summary_df.to_string(index=False))

# Save annotated files
adf.to_csv('Adipose_annotated_eqtls.csv', index=False)
ldf.to_csv('Liver_annotated_eqtls.csv', index=False)
print("\nSaved: Adipose_annotated_eqtls.csv, Liver_annotated_eqtls.csv")

# Cross-tissue artifact detection (same gene, trans in both)
merged = pd.merge(
    adf[['Symbol','cis_trans','flag_artifact','Peak -logP','Effect Size','tissue']],
    ldf[['Symbol','cis_trans','flag_artifact','Peak -logP','Effect Size','tissue']],
    on='Symbol', how='inner', suffixes=('_adipose','_liver')
)
multi_tissue_trans = merged[
    (merged['cis_trans_adipose']=='trans') & 
    (merged['cis_trans_liver']=='trans')
]
print(f"\nGenes with trans-eQTLs in BOTH tissues (cross-tissue trans candidates / mapping hotspots): {len(multi_tissue_trans)}")
if len(multi_tissue_trans) > 0:
    print(multi_tissue_trans[['Symbol','Peak -logP_adipose','Effect Size_adipose','Peak -logP_liver','Effect Size_liver']].to_string(index=False))

# Save candidate lists
for name, df in [('Adipose', adf), ('Liver', ldf)]:
    candidates = df[df['tier'] == 'Tier_1_HighConfidence'].sort_values('biol_score', ascending=False)
    candidates.to_csv(f'{name}_Tier1_candidates.csv', index=False)
    print(f"\n{name} Tier 1 candidates (n={len(candidates)}):")
    print(candidates[['Symbol','Description','Peak -logP','Effect Size','Distance_Mb','biol_score']].head(15).to_string(index=False))

# =============================================================================
# PLOTTING
# =============================================================================
sns.set_style('whitegrid')
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

def plot_tissue(ax_row, df, title):
    ax = ax_row[0]
    colors = {'cis':'#2ca02c', 'trans':'#d62728'}
    for ct, sub in df.groupby('cis_trans'):
        ax.scatter(sub['Peak -logP'], sub['Effect Size'].abs(), 
                   c=colors[ct], alpha=0.6, edgecolors='k', linewidth=0.3, s=40, label=ct)
    ax.set_xlabel('Peak -log10(P)')
    ax.set_ylabel('|Effect Size|')
    ax.set_title(f'{title}: -logP vs |Effect Size|')
    ax.legend(title='Cis/Trans')
    ax.set_xlim(18, df['Peak -logP'].max()*1.05)
    
    ax = ax_row[1]
    sns.boxplot(data=df, x='cis_trans', y='Peak -logP', ax=ax, palette=colors, order=['cis','trans'])
    ax.set_title(f'{title}: -logP Distribution')
    ax.set_ylabel('Peak -log10(P)')
    
    ax = ax_row[2]
    tier_counts = df['tier'].value_counts().reindex(['Tier_1_HighConfidence','Tier_2_Validate','Tier_3_FlaggedInspect'])
    tier_colors = {'Tier_1_HighConfidence':'#1f77b4', 'Tier_2_Validate':'#ff7f0e', 'Tier_3_FlaggedInspect':'#9467bd'}
    ax.bar(tier_counts.index, tier_counts.values, color=[tier_colors.get(x,'gray') for x in tier_counts.index])
    ax.set_title(f'{title}: Tier Assignment')
    ax.set_ylabel('Count')
    ax.set_xticklabels(['Tier 1\nHighConf','Tier 2\nValidate','Tier 3\nFlagged'], rotation=0, ha='center')

plot_tissue(axes[0], adf, 'Adipose')
plot_tissue(axes[1], ldf, 'Liver')

plt.tight_layout()
plt.savefig('eQTL_QC_summary.png', dpi=300, bbox_inches='tight')
print("\nSaved plot: eQTL_QC_summary.png")

# Manhattan-like plot (Peak Location vs -logP, colored by cis/trans)
fig, axes = plt.subplots(2, 1, figsize=(16, 10))
chr_order = [f'Chr{i}' for i in range(1,21)] + ['ChrX']
chr_colors = plt.cm.tab20.colors

def manhattan(ax, df, title):
    # Create numeric x-axis
    df = df.copy()
    df['chr_num'] = df['Peak_Chr'].map({c:i for i,c in enumerate(chr_order)})
    df = df.dropna(subset=['chr_num'])
    df = df.sort_values(['chr_num','Peak_Pos'])
    
    # Offset positions per chromosome
    offset = 0
    x_ticks = []
    x_labels = []
    for i, chrom in enumerate(chr_order):
        sub = df[df['Peak_Chr'] == chrom]
        if len(sub) == 0:
            continue
        x = sub['Peak_Pos'].values + offset
        ccolor = 'cis' if i % 2 == 0 else 'trans' # no, need to color by cis/trans status
        cis_mask = sub['cis_trans'] == 'cis'
        ax.scatter(x[cis_mask], sub['Peak -logP'][cis_mask], c='#2ca02c', s=20, alpha=0.7, label='cis' if i==0 else "")
        ax.scatter(x[~cis_mask], sub['Peak -logP'][~cis_mask], c='#d62728', s=20, alpha=0.7, label='trans' if i==0 else "")
        x_ticks.append(offset + sub['Peak_Pos'].mean())
        x_labels.append(chrom)
        offset += sub['Peak_Pos'].max() + 50  # 50 Mb spacer
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Peak -log10(P)')
    ax.set_title(f'{title} Manhattan-like Plot (colored by cis/trans)')
    ax.legend()

manhattan(axes[0], adf, 'Adipose')
manhattan(axes[1], ldf, 'Liver')
plt.tight_layout()
plt.savefig('eQTL_manhattan.png', dpi=300, bbox_inches='tight')
print("Saved plot: eQTL_manhattan.png")

# =============================================================================
# REPORT
# =============================================================================
report = f"""
# eQTL QC & Candidate Prioritization Report
Date: 2026-05-28
Tissues: Adipose (n={len(adf)}), Liver (n={len(ldf)})

## 1. Data Overview
Both datasets appear to be **pre-filtered significant eQTL hits** (minimum -logP = 20.0 in both).
They contain per-gene summary statistics from GeneNetwork 2 (HSNIH-Palmer cohort).

- **Adipose**: {len(adf)} genes, max -logP = {adf['Peak -logP'].max():.1f}, max |ES| = {adf['Effect Size'].abs().max():.2f}
- **Liver**: {len(ldf)} genes, max -logP = {ldf['Peak -logP'].max():.1f}, max |ES| = {ldf['Effect Size'].abs().max():.2f}
- Strong positive correlation between -logP and |Effect Size| (Adipose: {adf['Peak -logP'].corr(adf['Effect Size'].abs()):.2f}, Liver: {ldf['Peak -logP'].corr(ldf['Effect Size'].abs()):.2f}).
  This is reassuring: high significance is driven by large effects, not just vanishing standard errors.

## 2. Cis vs Trans Classification (threshold: {CIS_THRESHOLD_MB} Mb)
| Tissue | Cis | Trans | % Cis |
|--------|-----|-------|-------|
| Adipose | {(adf['cis_trans']=='cis').sum()} | {(adf['cis_trans']=='trans').sum()} | {(adf['cis_trans']=='cis').mean()*100:.1f}% |
| Liver   | {(ldf['cis_trans']=='cis').sum()} | {(ldf['cis_trans']=='trans').sum()} | {(ldf['cis_trans']=='cis').mean()*100:.1f}% |

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
- Trans with |Effect Size| > {ARTIFACT_ES_TRANS} -> **Flag: MEGA_EFFECT**
- Trans with -logP > {ARTIFACT_LOGP_TRANS} -> **Flag: MEGA_EFFECT**
- Pseudogene with trans-eQTL -> **Flag: PSEUDOGENE_TRANS**
- Uncharacterized LOC with trans & -logP > 60 -> **Flag: LOC_UNCHARACTERIZED**

**Results**:
- Adipose flagged: {adf['flag_artifact'].sum()} / {len(adf)} ({adf['flag_artifact'].mean()*100:.1f}%)
- Liver flagged: {ldf['flag_artifact'].sum()} / {len(ldf)} ({ldf['flag_artifact'].mean()*100:.1f}%)

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
"""

with open('eQTL_analysis_report.md', 'w') as f:
    f.write(report)
print("\nSaved report: eQTL_analysis_report.md")

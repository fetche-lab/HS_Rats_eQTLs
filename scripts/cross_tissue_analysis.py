import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
# from matplotlib_venn import venn2  # not installed, using alternative
import warnings
warnings.filterwarnings('ignore')

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300

# =============================================================================
# LOAD
# =============================================================================
adf = pd.read_csv('Adipose_annotated_eqtls.csv')
ldf = pd.read_csv('Liver_annotated_eqtls.csv')

# =============================================================================
# REFINED GENOMIC CONTEXT
# =============================================================================
def genomic_context(row):
    if row['cis_trans'] == 'trans':
        return 'trans'
    d = row['Distance_Mb']
    if d < 0.1:
        return 'promoter'
    elif d < 1.0:
        return 'proximal'
    else:
        return 'distal_cis'

for df in [adf, ldf]:
    df['genomic_context'] = df.apply(genomic_context, axis=1)

# =============================================================================
# CROSS-TISSUE MERGE
# =============================================================================
# Use outer merge so we capture tissue-specific genes too
merged = pd.merge(
    adf, ldf,
    on='Symbol', how='outer', suffixes=('_adipose', '_liver'),
    indicator=True
)
merged['in_adipose'] = merged['_merge'].isin(['left_only', 'both'])
merged['in_liver'] = merged['_merge'].isin(['right_only', 'both'])
merged['in_both'] = merged['_merge'] == 'both'

# =============================================================================
# SHARED PEAK & DIRECTION ANALYSIS (for genes in both)
# =============================================================================
merged['shared_peak'] = False
merged['same_direction'] = False
merged['peak_distance_mb'] = np.nan

both_mask = merged['in_both']
if both_mask.sum() > 0:
    same_chr = merged.loc[both_mask, 'Peak_Chr_adipose'] == merged.loc[both_mask, 'Peak_Chr_liver']
    peak_dist = np.abs(merged.loc[both_mask, 'Peak_Pos_adipose'].fillna(0) - merged.loc[both_mask, 'Peak_Pos_liver'].fillna(0))
    merged.loc[both_mask, 'peak_distance_mb'] = peak_dist
    merged.loc[both_mask, 'shared_peak'] = same_chr.values & (peak_dist.values <= 4.0)
    
    es_prod = merged.loc[both_mask, 'Effect Size_adipose'].fillna(0) * merged.loc[both_mask, 'Effect Size_liver'].fillna(0)
    merged.loc[both_mask, 'same_direction'] = es_prod.values > 0

# =============================================================================
# PLEIOTROPY CLASSIFICATION
# =============================================================================
def classify_pleiotropy(row):
    if not row['in_both']:
        return 'tissue_specific'
    
    cis_a = row['cis_trans_adipose'] == 'cis' if pd.notna(row['cis_trans_adipose']) else False
    cis_l = row['cis_trans_liver'] == 'cis' if pd.notna(row['cis_trans_liver']) else False
    trans_a = row['cis_trans_adipose'] == 'trans' if pd.notna(row['cis_trans_adipose']) else False
    trans_l = row['cis_trans_liver'] == 'trans' if pd.notna(row['cis_trans_liver']) else False
    
    if trans_a and trans_l:
        return 'trans_both_hotspot'  # renamed: not assumed artifact, but flagged for inspection
    
    if cis_a and cis_l and row['shared_peak']:
        if row['same_direction']:
            return 'shared_cis_conserved'
        else:
            return 'shared_cis_flip'
    
    if cis_a and cis_l and not row['shared_peak']:
        return 'shared_cis_different_locus'
    
    if cis_a or cis_l:
        return 'tissue_specific_cis'
    
    return 'other'

merged['pleiotropy_class'] = merged.apply(classify_pleiotropy, axis=1)

# =============================================================================
# TISSUE SPECIFICITY SCORE (for genes in both)
# =============================================================================
# Compute a specificity index based on normalized biological scores
# 0 = equally strong in both tissues, 1 = completely specific to one tissue
both = merged[merged['in_both']].copy()
if len(both) > 0:
    ba = both['biol_score_adipose'].fillna(0)
    bl = both['biol_score_liver'].fillna(0)
    merged.loc[merged['in_both'], 'tissue_specificity_index'] = np.abs(ba - bl) / (ba + bl + 0.01)
    
    # Also compute a "conservation score" = average biol score * (1 - specificity/2)
    # High when strong in both tissues, low when specific or weak
    avg_score = (ba + bl) / 2
    conservation = avg_score * (1 - merged.loc[merged['in_both'], 'tissue_specificity_index'] / 2)
    merged.loc[merged['in_both'], 'conservation_score'] = conservation
else:
    merged['tissue_specificity_index'] = np.nan
    merged['conservation_score'] = np.nan

# =============================================================================
# MASTER CANDIDATE SCORING
# =============================================================================
def master_score(row):
    """
    Score candidates for follow-up.
    Higher = better cross-tissue evidence or very strong tissue-specific signal.
    """
    score = 0.0
    
    # Base: best biol score across tissues
    biol_scores = []
    if pd.notna(row['biol_score_adipose']): biol_scores.append(row['biol_score_adipose'])
    if pd.notna(row['biol_score_liver']): biol_scores.append(row['biol_score_liver'])
    score += max(biol_scores) * 40  # 0-40 pts
    
    # Cross-tissue bonuses
    if row['in_both']:
        if row['pleiotropy_class'] == 'shared_cis_conserved':
            score += 30  # Strong pleiotropic conserved signal
            if row['same_direction'] and pd.notna(row['Effect Size_adipose']) and pd.notna(row['Effect Size_liver']):
                # Bonus if both effect sizes are large
                min_es = min(abs(row['Effect Size_adipose']), abs(row['Effect Size_liver']))
                score += min_es * 5
        elif row['pleiotropy_class'] == 'shared_cis_flip':
            score += 15  # Interesting but needs explanation
        elif row['pleiotropy_class'] == 'shared_cis_different_locus':
            score += 10  # Different regulatory mechanisms
        elif row['pleiotropy_class'] == 'tissue_specific_cis':
            score += 12  # Tissue-specific regulation is interesting
    else:
        score += 5  # Tissue-specific gene
    
    # Penalties
    if row['pleiotropy_class'] == 'trans_both_hotspot':
        score -= 20  # reduced penalty: these are kept for inspection, not discarded
    
    # Bonus for promoter-proximal peaks (more interpretable)
    contexts = []
    if pd.notna(row.get('genomic_context_adipose')): contexts.append(row['genomic_context_adipose'])
    if pd.notna(row.get('genomic_context_liver')): contexts.append(row['genomic_context_liver'])
    if 'promoter' in contexts:
        score += 5
    elif 'proximal' in contexts:
        score += 3
    
    return score

merged['master_score'] = merged.apply(master_score, axis=1)

# =============================================================================
# TRANS-BAND CLUSTERING (new) — MOVED BEFORE REPORT GENERATION
# =============================================================================
print("\n--- Clustering trans-bands ---")
trans_all = pd.concat([
    adf[adf['cis_trans'] == 'trans'][['Symbol', 'Peak_Chr', 'Peak_Pos', 'Peak -logP', 'Effect Size', 'trans_flag']].assign(Tissue='Adipose'),
    ldf[ldf['cis_trans'] == 'trans'][['Symbol', 'Peak_Chr', 'Peak_Pos', 'Peak -logP', 'Effect Size', 'trans_flag']].assign(Tissue='Liver')
])
trans_all['Band_Mb'] = (trans_all['Peak_Pos'] / 5).astype(int) * 5
trans_all['Band_Label'] = trans_all['Peak_Chr'].astype(str) + ':' + trans_all['Band_Mb'].astype(str) + '-' + (trans_all['Band_Mb'] + 5).astype(str) + ' Mb'

# Reset index to avoid duplicate loc issues
trans_all = trans_all.reset_index(drop=True)

band_summary = trans_all.groupby(['Peak_Chr', 'Band_Mb']).agg(
    Band_Label=('Band_Label', 'first'),
    N_Genes=('Symbol', 'nunique'),
    Tissues=('Tissue', lambda x: ', '.join(sorted(x.unique()))),
    Top_Gene=('Peak -logP', lambda x: trans_all.loc[x.idxmax(), 'Symbol']),
    Max_logP=('Peak -logP', 'max'),
    Max_ES=('Effect Size', lambda x: trans_all.loc[x.abs().idxmax(), 'Effect Size']),
    Genes=('Symbol', lambda x: ', '.join(sorted(x.unique())))
).reset_index()

band_summary['Cross_Tissue'] = band_summary['Tissues'].str.contains('Adipose') & band_summary['Tissues'].str.contains('Liver')
band_summary['Priority'] = band_summary['Cross_Tissue'].astype(int) * 1000 + band_summary['N_Genes'] * 10 + band_summary['Max_logP']
band_summary = band_summary.sort_values('Priority', ascending=False)

def annotate_band(row):
    flags = []
    genes = row['Genes'].split(', ')
    if row['Cross_Tissue']:
        flags.append('CROSS_TISSUE_TRANS')
    if row['Max_logP'] > 80 or abs(float(row['Max_ES'])) > 2.0:
        flags.append('MEGA_EFFECT')
    if any('-ps' in g for g in genes):
        flags.append('PSEUDOGENE_TRANS')
    if any(g.startswith('LOC') or g.startswith('RGD') for g in genes):
        flags.append('LOC_UNCHARACTERIZED')
    if row['N_Genes'] >= 3:
        flags.append('MASTER_REGULATOR_CANDIDATE')
    if row['Band_Mb'] < 2:
        flags.append('NEAR_TELOMERE')
    if any('Rnase' in g for g in genes):
        flags.append('PARALOG_FAMILY')
    return '; '.join(flags)

band_summary['Flags'] = band_summary.apply(annotate_band, axis=1)
band_out = band_summary[['Band_Label', 'N_Genes', 'Tissues', 'Top_Gene', 'Max_logP', 'Max_ES', 'Cross_Tissue', 'Flags', 'Genes']]
band_out.to_csv('trans_bands.csv', index=False)
print(f"Saved trans_bands.csv with {len(band_out)} bands")

# Trans-band landscape plot
fig, ax = plt.subplots(figsize=(14, 6))
chr_order = [f'Chr{i}' for i in range(1, 21)] + ['ChrX']
all_data = pd.concat([adf, ldf])
chr_offsets = {}
offset = 0
for c in chr_order:
    max_pos = all_data[all_data['Peak_Chr'] == c]['Peak_Pos'].max()
    if pd.isna(max_pos):
        max_pos = 200
    chr_offsets[c] = offset
    offset += max_pos + 10

colors = {'Adipose': '#E69F00', 'Liver': '#56B4E9', 'Adipose, Liver': '#CC79A7'}
for _, row in band_out.iterrows():
    c = row['Band_Label'].split(':')[0]
    mb = int(row['Band_Label'].split(':')[1].split('-')[0])
    x = chr_offsets.get(c, 0) + mb
    y = row['N_Genes']
    size = row['Max_logP'] * 3
    color = colors.get(row['Tissues'], '#999999')
    alpha = 0.9 if row['Cross_Tissue'] else 0.5
    ax.scatter(x, y, s=size, c=color, alpha=alpha, edgecolors='black', linewidth=0.5)
    if row['N_Genes'] >= 3 or row['Cross_Tissue']:
        ax.annotate(row['Band_Label'].replace('Chr', ''), (x, y), textcoords='offset points', xytext=(0, 8), ha='center', fontsize=7)

for c in chr_order:
    if c in chr_offsets:
        ax.axvline(chr_offsets[c] - 10, color='gray', linestyle='--', alpha=0.3)
        ax.text(chr_offsets[c] + (all_data[all_data['Peak_Chr'] == c]['Peak_Pos'].max() or 0)/2, -0.8, c.replace('Chr', ''), ha='center', fontsize=9)

ax.set_ylabel('Number of trans-eQTL genes in 5 Mb band')
ax.set_xlabel('Chromosome')
ax.set_title('Trans-Band Landscape: Genomic Hotspots of Distal Regulation\n(Point size = max -logP; Purple = cross-tissue; Orange = adipose; Blue = liver)')
ax.set_xlim(-20, offset)
ax.set_ylim(-1, band_out['N_Genes'].max() + 2)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E69F00', markersize=10, label='Adipose only'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#56B4E9', markersize=10, label='Liver only'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#CC79A7', markersize=10, label='Cross-tissue'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=15, label='Large -logP')
]
ax.legend(handles=legend_elements, loc='upper right')
plt.tight_layout()
plt.savefig('trans_band_landscape.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: trans_band_landscape.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
summary_lines = []
summary_lines.append("# Cross-Tissue eQTL Integration Report")
summary_lines.append(f"Date: 2026-05-28")
summary_lines.append("")
summary_lines.append("## 1. Dataset Overlap")
summary_lines.append(f"- Adipose-only genes: {(merged['in_adipose'] & ~merged['in_liver']).sum()}")
summary_lines.append(f"- Liver-only genes: {(merged['in_liver'] & ~merged['in_adipose']).sum()}")
summary_lines.append(f"- Genes in both tissues: {merged['in_both'].sum()}")
summary_lines.append("")

summary_lines.append("## 2. Genomic Context Distribution")
for tissue, df in [('Adipose', adf), ('Liver', ldf)]:
    summary_lines.append(f"\n### {tissue}")
    ctx = df['genomic_context'].value_counts()
    for ctx_name, count in ctx.items():
        pct = count / len(df) * 100
        summary_lines.append(f"- {ctx_name}: {count} ({pct:.1f}%)")

summary_lines.append("")
summary_lines.append("## 3. Pleiotropy Classification (genes in both tissues)")
pleio_counts = merged[merged['in_both']]['pleiotropy_class'].value_counts()
for cls, count in pleio_counts.items():
    summary_lines.append(f"- {cls}: {count}")

summary_lines.append("")
summary_lines.append("## 3b. Trans-Band Hotspots")
summary_lines.append(f"Total trans-bands identified: {len(band_out)}")
for _, row in band_out.head(10).iterrows():
    ct = " [CROSS-TISSUE]" if row['Cross_Tissue'] else ""
    summary_lines.append(f"- {row['Band_Label']}: {row['N_Genes']} genes, top {row['Top_Gene']} (-logP={row['Max_logP']:.1f}){ct}")
    summary_lines.append(f"  Flags: {row['Flags']}")

summary_lines.append("")
summary_lines.append("## 4. Shared Conserved cis-eQTLs (same peak, same direction)")
conserved = merged[merged['pleiotropy_class'] == 'shared_cis_conserved'].sort_values('master_score', ascending=False)
summary_lines.append(f"Count: {len(conserved)}")
if len(conserved) > 0:
    summary_lines.append("")
    summary_lines.append("| Symbol | Description | Adipose -logP | Adipose ES | Liver -logP | Liver ES | Peak Dist (Mb) |")
    summary_lines.append("|--------|-------------|---------------|------------|-------------|----------|----------------|")
    for _, row in conserved.head(20).iterrows():
        desc = str(row['Description_adipose']) if pd.notna(row['Description_adipose']) else str(row['Description_liver'])
        desc = desc[:50]
        summary_lines.append(f"| {row['Symbol']} | {desc} | {row['Peak -logP_adipose']:.1f} | {row['Effect Size_adipose']:.2f} | {row['Peak -logP_liver']:.1f} | {row['Effect Size_liver']:.2f} | {row['peak_distance_mb']:.3f} |")

summary_lines.append("")
summary_lines.append("## 5. Tissue-Specific Strong cis-eQTLs")
ts_cis = merged[merged['pleiotropy_class'] == 'tissue_specific_cis'].sort_values('master_score', ascending=False)
summary_lines.append(f"Count: {len(ts_cis)}")
if len(ts_cis) > 0:
    summary_lines.append("")
    summary_lines.append("| Symbol | Tissue | -logP | ES | Context | Description |")
    summary_lines.append("|--------|--------|-------|----|---------|-------------|")
    for _, row in ts_cis.head(15).iterrows():
        tissue = 'Adipose' if row['cis_trans_adipose']=='cis' else 'Liver'
        lp = row['Peak -logP_adipose'] if tissue=='Adipose' else row['Peak -logP_liver']
        es = row['Effect Size_adipose'] if tissue=='Adipose' else row['Effect Size_liver']
        ctx = row['genomic_context_adipose'] if tissue=='Adipose' else row['genomic_context_liver']
        desc = str(row['Description_adipose']) if tissue=='Adipose' else str(row['Description_liver'])
        desc = desc[:45]
        summary_lines.append(f"| {row['Symbol']} | {tissue} | {lp:.1f} | {es:.2f} | {ctx} | {desc} |")

summary_lines.append("")
summary_lines.append("## 6. Top Master Candidates (cross-tissue scored)")
top_master = merged[merged['pleiotropy_class'] != 'trans_both_hotspot'].sort_values('master_score', ascending=False)
summary_lines.append(f"Total scored: {len(top_master)}")
summary_lines.append("")
summary_lines.append("| Rank | Symbol | Pleiotropy | Adipose Tier | Liver Tier | Master Score |")
summary_lines.append("|------|--------|------------|--------------|------------|--------------|")
for i, (_, row) in enumerate(top_master.head(25).iterrows(), 1):
    a_tier = row['tier_adipose'] if pd.notna(row.get('tier_adipose')) else 'N/A'
    l_tier = row['tier_liver'] if pd.notna(row.get('tier_liver')) else 'N/A'
    summary_lines.append(f"| {i} | {row['Symbol']} | {row['pleiotropy_class']} | {a_tier} | {l_tier} | {row['master_score']:.1f} |")

summary_text = "\n".join(summary_lines)
with open('cross_tissue_report.md', 'w') as f:
    f.write(summary_text)
print("Saved: cross_tissue_report.md")

# =============================================================================
# SAVE TABLES
# =============================================================================
merged.to_csv('cross_tissue_merged.csv', index=False)
print("Saved: cross_tissue_merged.csv")

top_master.to_csv('cross_tissue_master_candidates.csv', index=False)
print("Saved: cross_tissue_master_candidates.csv")

# =============================================================================
# PLOTTING
# =============================================================================

# --- Plot 1: Effect Size Comparison ---
fig, ax = plt.subplots(figsize=(8, 8))
both_df = merged[merged['in_both']].copy()
# NOW INCLUDING trans_both_hotspot, flagged with distinct markers
plot_df = both_df.copy()

pleio_colors = {
    'shared_cis_conserved': '#1f77b4',
    'shared_cis_flip': '#ff7f0e',
    'shared_cis_different_locus': '#2ca02c',
    'tissue_specific_cis': '#d62728',
    'trans_both_hotspot': '#000000',
    'other': '#9467bd'
}
pleio_markers = {
    'shared_cis_conserved': 'o',
    'shared_cis_flip': 'o',
    'shared_cis_different_locus': 'o',
    'tissue_specific_cis': 'o',
    'trans_both_hotspot': '^',
    'other': 'o'
}

for cls, sub in plot_df.groupby('pleiotropy_class'):
    ax.scatter(sub['Effect Size_adipose'], sub['Effect Size_liver'],
               c=pleio_colors.get(cls, 'gray'), marker=pleio_markers.get(cls, 'o'),
               label=cls, alpha=0.75, edgecolors='k', linewidth=0.3, s=70)

ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
ax.axvline(0, color='black', linestyle='--', linewidth=0.5)
ax.plot([-5, 5], [-5, 5], 'k--', linewidth=0.5, alpha=0.5, label='y=x')
ax.set_xlabel('Adipose Effect Size')
ax.set_ylabel('Liver Effect Size')
ax.set_title('Cross-Tissue Effect Size Comparison\n(genes detected in both tissues; black triangles = trans hotspots)')
ax.legend(loc='upper left', fontsize=8)
ax.set_xlim(-5.5, 5.5)
ax.set_ylim(-5.5, 5.5)
plt.tight_layout()
plt.savefig('cross_tissue_effectsize_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cross_tissue_effectsize_comparison.png")

# --- Plot 2: -logP Comparison ---
fig, ax = plt.subplots(figsize=(8, 8))
for cls, sub in plot_df.groupby('pleiotropy_class'):
    ax.scatter(sub['Peak -logP_adipose'], sub['Peak -logP_liver'],
               c=pleio_colors.get(cls, 'gray'), marker=pleio_markers.get(cls, 'o'),
               label=cls, alpha=0.75, edgecolors='k', linewidth=0.3, s=70)

ax.plot([20, 180], [20, 180], 'k--', linewidth=0.5, alpha=0.5)
ax.set_xlabel('Adipose Peak -log10(P)')
ax.set_ylabel('Liver Peak -log10(P)')
ax.set_title('Cross-Tissue Significance Comparison\n(black triangles = trans hotspots)')
ax.legend(loc='upper left', fontsize=8)
ax.set_xlim(18, 180)
ax.set_ylim(18, 180)
plt.tight_layout()
plt.savefig('cross_tissue_logp_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cross_tissue_logp_comparison.png")

# --- Plot 3: Pleiotropy Class Bar Chart ---
fig, ax = plt.subplots(figsize=(10, 5))
pleio_counts = plot_df['pleiotropy_class'].value_counts()
colors = [pleio_colors.get(k, 'gray') for k in pleio_counts.index]
ax.barh(pleio_counts.index, pleio_counts.values, color=colors, edgecolor='black')
ax.set_xlabel('Number of Genes')
ax.set_title('Pleiotropy Classification (genes in both tissues)\nblack = trans_both_hotspot')
for i, v in enumerate(pleio_counts.values):
    ax.text(v + 1, i, str(v), va='center')
plt.tight_layout()
plt.savefig('cross_tissue_pleiotropy.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cross_tissue_pleiotropy.png")

# --- Plot 4: Genomic Context ---
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
ctx_order = ['promoter', 'proximal', 'distal_cis', 'trans']
ctx_palette = {'promoter':'#e41a1c', 'proximal':'#377eb8', 'distal_cis':'#4daf4a', 'trans':'#984ea3'}

for ax, (tissue, df) in zip(axes, [('Adipose', adf), ('Liver', ldf)]):
    ctx_counts = df['genomic_context'].value_counts().reindex(ctx_order, fill_value=0)
    ax.bar(ctx_counts.index, ctx_counts.values, color=[ctx_palette.get(x,'gray') for x in ctx_counts.index], edgecolor='black')
    ax.set_title(f'{tissue} Genomic Context')
    ax.set_ylabel('Count')
    ax.set_xticklabels(ctx_counts.index, rotation=30, ha='right')
    for i, v in enumerate(ctx_counts.values):
        ax.text(i, v + 5, str(v), ha='center', va='bottom')

plt.tight_layout()
plt.savefig('cross_tissue_genomic_context.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cross_tissue_genomic_context.png")

# --- Plot 5: Tissue Specificity Index ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Histogram
ax = axes[0]
spec_vals = merged['tissue_specificity_index'].dropna()
ax.hist(spec_vals, bins=20, color='steelblue', edgecolor='black', alpha=0.8)
ax.set_xlabel('Tissue Specificity Index')
ax.set_ylabel('Number of Genes')
ax.set_title('Distribution of Tissue Specificity\n(0 = equal strength, 1 = tissue-specific)')
ax.axvline(spec_vals.median(), color='red', linestyle='--', label=f'Median = {spec_vals.median():.2f}')
ax.legend()

# Boxplot by pleiotropy class
ax = axes[1]
box_data = []
box_labels = []
for cls in ['shared_cis_conserved', 'shared_cis_flip', 'tissue_specific_cis', 'shared_cis_different_locus', 'trans_both_hotspot']:
    vals = merged[merged['pleiotropy_class'] == cls]['tissue_specificity_index'].dropna()
    if len(vals) > 0:
        box_data.append(vals)
        box_labels.append(cls.replace('_', '\n'))

bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)
for patch, color in zip(bp['boxes'], [pleio_colors.get(l.replace('\n','_'), 'gray') for l in box_labels]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_ylabel('Tissue Specificity Index')
ax.set_title('Tissue Specificity by Pleiotropy Class')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right')

plt.tight_layout()
plt.savefig('cross_tissue_specificity.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cross_tissue_specificity.png")

# --- Plot 6: Tier 1 Overlap Bar Chart (alternative to venn) ---
fig, ax = plt.subplots(figsize=(8, 6))
t1_adipose = set(adf[adf['tier'] == 'Tier_1_HighConfidence']['Symbol'])
t1_liver = set(ldf[ldf['tier'] == 'Tier_1_HighConfidence']['Symbol'])
only_a = len(t1_adipose - t1_liver)
only_l = len(t1_liver - t1_adipose)
both_t1 = len(t1_adipose & t1_liver)

categories = ['Adipose only', 'Both', 'Liver only']
values = [only_a, both_t1, only_l]
colors = ['#ff7f0e', '#9467bd', '#2ca02c']
bars = ax.bar(categories, values, color=colors, edgecolor='black')
ax.set_ylabel('Number of Genes')
ax.set_title('Tier 1 Candidate Overlap')
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, str(val), ha='center', va='bottom', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('cross_tissue_tier1_overlap.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cross_tissue_tier1_overlap.png")

# =============================================================================
# FINAL CONSOLE SUMMARY
# =============================================================================
print("\n" + "="*60)
print("CROSS-TISSUE ANALYSIS COMPLETE")
print("="*60)
print(f"Genes in both tissues: {merged['in_both'].sum()}")
print(f"Shared conserved cis-eQTLs: {len(conserved)}")
print(f"Tissue-specific cis-eQTLs: {len(ts_cis)}")
print(f"Cross-tissue trans hotspots: {len(merged[merged['pleiotropy_class']=='trans_both_hotspot'])}")
print(f"Trans-bands identified: {len(band_out)}")
print(f"\nTop 10 Master Candidates:")
print(top_master[['Symbol','pleiotropy_class','master_score']].head(10).to_string(index=False))

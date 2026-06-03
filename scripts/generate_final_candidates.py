import pandas as pd
import numpy as np

# =============================================================================
# LOAD CROSS-TISSUE DATA
# =============================================================================
merged = pd.read_csv('cross_tissue_merged.csv')

# =============================================================================
# FILTER TO KNOWN GENES ONLY
# =============================================================================
known_mask = (
    ~merged['Symbol'].str.startswith('LOC', na=False) &
    ~merged['Symbol'].str.contains('-ps', na=False) &  # pseudogenes
    ~merged['Description_adipose'].fillna('').str.contains('uncharacterized', case=False) &
    ~merged['Description_liver'].fillna('').str.contains('uncharacterized', case=False) &
    ~merged['Description_adipose'].fillna('').str.contains('pseudogene', case=False) &
    ~merged['Description_liver'].fillna('').str.contains('pseudogene', case=False)
)
known = merged[known_mask].copy()

print(f"Total merged genes: {len(merged)}")
print(f"Known genes: {len(known)}")

# =============================================================================
# FUNCTIONAL CATEGORY ANNOTATION
# =============================================================================
def get_description(row):
    """Pick the best available description."""
    if pd.notna(row['Description_adipose']) and row['Description_adipose'] != '':
        return str(row['Description_adipose'])
    if pd.notna(row['Description_liver']) and row['Description_liver'] != '':
        return str(row['Description_liver'])
    return ''

known['best_description'] = known.apply(get_description, axis=1)

def classify_function(desc):
    desc = desc.lower()
    scores = {}
    
    # Lipid/Metabolic
    scores['lipid_metabolism'] = sum([
        'lipid' in desc, 'fatty acid' in desc, 'acyl' in desc, 'coa' in desc,
        'desaturase' in desc, 'synthase' in desc, 'sterol' in desc,
        'glucose' in desc, 'insulin' in desc, 'bile' in desc,
        ' UDP' in desc, 'glucuronosyltransferase' in desc, 'sulfotransferase' in desc,
        'cytochrome P450' in desc, 'glutathione' in desc, 'peroxisome' in desc,
        'mitochondria' in desc, 'oxidoreductase' in desc, 'dehydrogenase' in desc,
        'lipase' in desc, 'phospholipase' in desc, 'transferase' in desc,
        'acyltransferase' in desc, 'amino acid' in desc, 'carbohydrate' in desc
    ])
    
    # Immune / Inflammation
    scores['immune'] = sum([
        'interferon' in desc, 'interleukin' in desc, 'cytokine' in desc,
        'chemokine' in desc, 'tumor necrosis factor' in desc, 'TNF' in desc,
        'MHC' in desc, 'histocompatibility' in desc, 'RT1 class' in desc,
        'immunoglobulin' in desc, 'lymphocyte' in desc, 'leukocyte' in desc,
        'antigen' in desc, 'complement' in desc, 'defensin' in desc,
        ' major histocompatibility' in desc
    ])
    
    # Transcription / Chromatin
    scores['transcription_chromatin'] = sum([
        'transcription' in desc, 'transcriptional' in desc,
        'histone' in desc, 'chromatin' in desc, 'zinc finger' in desc,
        'homeobox' in desc, 'helix-loop-helix' in desc, 'bHLH' in desc,
        'nuclear receptor' in desc, 'TFII' in desc
    ])
    
    # Signal transduction
    scores['signaling'] = sum([
        'receptor' in desc, 'kinase' in desc, 'phosphatase' in desc,
        'G protein' in desc, 'GTPase' in desc, 'phosphodiesterase' in desc,
        'cyclase' in desc, 'phosphoinositide' in desc
    ])
    
    # Proteolysis / Protein turnover
    scores['proteolysis'] = sum([
        'protease' in desc, 'peptidase' in desc, 'ubiquitin' in desc,
        'proteasome' in desc, 'caspase' in desc, 'cathepsin' in desc
    ])
    
    # Extracellular / Matrix
    scores['extracellular_matrix'] = sum([
        'collagen' in desc, 'laminin' in desc, 'fibronectin' in desc,
        'matrix metalloproteinase' in desc, 'gelatinase' in desc
    ])
    
    best = max(scores, key=scores.get)
    if scores[best] == 0:
        return 'other'
    return best

known['functional_category'] = known['best_description'].apply(classify_function)

# =============================================================================
# BUILD FINAL ACTIONABLE TABLE
# =============================================================================
# Focus on high-confidence: either Tier 1 in at least one tissue, or shared_cis_conserved with good scores
high_conf = known[
    (known['tier_adipose'] == 'Tier_1_HighConfidence') |
    (known['tier_liver'] == 'Tier_1_HighConfidence') |
    (
        (known['pleiotropy_class'] == 'shared_cis_conserved') &
        (known['master_score'] > 50)
    )
].copy()

print(f"High-confidence known candidates: {len(high_conf)}")

# =============================================================================
# PANGENOMIC CANDIDATES (NEW): flagged trans + uncharacterized genes
# =============================================================================
pangeno_mask = (
    (merged['pleiotropy_class'] == 'trans_both_hotspot') |
    (merged['Symbol'].str.startswith('LOC', na=False)) |
    (merged['Symbol'].str.contains('-ps', na=False)) |
    (merged['Description_adipose'].fillna('').str.contains('uncharacterized', case=False)) |
    (merged['Description_liver'].fillna('').str.contains('uncharacterized', case=False))
)
pangeno = merged[pangeno_mask].copy()
print(f"Pangenomic candidates (flagged trans / uncharacterized): {len(pangeno)}")
if len(pangeno) > 0:
    pcols = ['Symbol', 'pleiotropy_class', 'master_score', 'Peak -logP_adipose', 'Peak -logP_liver',
             'Effect Size_adipose', 'Effect Size_liver', 'Peak_Chr_adipose', 'Peak_Chr_liver',
             'Peak_Pos_adipose', 'Peak_Pos_liver']
    pcols = [c for c in pcols if c in pangeno.columns]
    pangeno[pangeno.columns[pangeno.columns.str.contains('trans_flag')]].head()
    pangeno[pangeno.columns[pangeno.columns.str.contains('tier')]].head()
    pangeno_out = pangeno[['Symbol'] + [c for c in pcols if c != 'Symbol']].sort_values('master_score', ascending=False)
    pangeno_out.to_csv('pangenomic_candidates.csv', index=False)
    print("Saved: pangenomic_candidates.csv")

# Select columns for the final clean output
cols = [
    'Symbol', 'best_description', 'functional_category',
    'pleiotropy_class', 'master_score',
    'Peak -logP_adipose', 'Effect Size_adipose', 'genomic_context_adipose', 'tier_adipose',
    'Peak -logP_liver', 'Effect Size_liver', 'genomic_context_liver', 'tier_liver',
    'same_direction', 'peak_distance_mb'
]

# Only keep columns that exist
final_cols = [c for c in cols if c in high_conf.columns]
final_df = high_conf[final_cols].sort_values('master_score', ascending=False)

# Round floats for readability
for col in final_df.columns:
    if final_df[col].dtype == 'float64':
        final_df[col] = final_df[col].round(3)

final_df.to_csv('final_actionable_candidates.csv', index=False)
print("Saved: final_actionable_candidates.csv")

# =============================================================================
# SUMMARY BY CATEGORY
# =============================================================================
cat_summary = final_df['functional_category'].value_counts()
print("\n=== Functional Category Breakdown ===")
print(cat_summary.to_string())

# =============================================================================
# WRITE FINAL RECOMMENDATIONS REPORT
# =============================================================================
report = []
report.append("# Final Actionable Candidate Recommendations")
report.append("")
report.append("## Filtering Applied")
report.append("- Excluded: LOC genes, uncharacterized transcripts, pseudogenes")
report.append("- Included: Tier 1 in at least one tissue OR shared conserved cis-eQTL with master score > 50")
report.append(f"- Result: **{len(final_df)}** high-confidence known genes for follow-up")
report.append("")
report.append("## Top Candidates by Biological Theme")
report.append("")

for cat in cat_summary.index:
    report.append(f"### {cat.replace('_', ' ').title()}")
    sub = final_df[final_df['functional_category'] == cat].head(10)
    report.append("| Symbol | Description | Pleiotropy | Adipose -logP | Liver -logP | Master Score |")
    report.append("|--------|-------------|------------|---------------|-------------|--------------|")
    for _, row in sub.iterrows():
        desc = row['best_description'][:55]
        a_lp = f"{row['Peak -logP_adipose']:.1f}" if pd.notna(row['Peak -logP_adipose']) else 'N/A'
        l_lp = f"{row['Peak -logP_liver']:.1f}" if pd.notna(row['Peak -logP_liver']) else 'N/A'
        report.append(f"| {row['Symbol']} | {desc} | {row['pleiotropy_class']} | {a_lp} | {l_lp} | {row['master_score']:.1f} |")
    report.append("")

report.append("## Pangenomic Candidates (Flagged Trans / Uncharacterized)")
report.append(f"A separate list of **{len(pangeno)}** genes with trans-eQTLs in both tissues, pseudogenes, or uncharacterized LOC symbols is saved as `pangenomic_candidates.csv`.")
report.append("These are **not discarded** but are prioritized for:")
report.append("- Pangenomic graph re-mapping (e.g., Minigraph-Cactus with 8 HS founder genomes)")
report.append("- Structural variant / CNV inspection in RGD")
report.append("- Read-depth analysis from raw RNA-seq BAMs")
report.append("- Founder haplotype tracing to confirm allele-specific segregation")
report.append("")

report.append("## Priority Tiers for Experimental Follow-up")
report.append("")
report.append("### Immediate Priority (shared conserved cis, Tier 1 in both tissues)")
top_shared = final_df[
    (final_df['pleiotropy_class'] == 'shared_cis_conserved') &
    (final_df['tier_adipose'] == 'Tier_1_HighConfidence') &
    (final_df['tier_liver'] == 'Tier_1_HighConfidence')
].head(15)
for _, row in top_shared.iterrows():
    report.append(f"- **{row['Symbol']}** ({row['functional_category'].replace('_',' ')}): master score {row['master_score']:.1f}")
report.append("")

report.append("### High Priority (Tier 1 in one tissue, strong in other)")
top_one = final_df[
    (
        ((final_df['tier_adipose'] == 'Tier_1_HighConfidence') & (final_df['tier_liver'] != 'Tier_1_HighConfidence')) |
        ((final_df['tier_liver'] == 'Tier_1_HighConfidence') & (final_df['tier_adipose'] != 'Tier_1_HighConfidence'))
    ) &
    (final_df['pleiotropy_class'] == 'shared_cis_conserved')
].head(15)
for _, row in top_one.iterrows():
    report.append(f"- **{row['Symbol']}** ({row['functional_category'].replace('_',' ')}): master score {row['master_score']:.1f}")
report.append("")

report.append("### Validate (shared conserved but Tier 2 in both)")
top_t2 = final_df[
    (final_df['pleiotropy_class'] == 'shared_cis_conserved') &
    (final_df['tier_adipose'] != 'Tier_1_HighConfidence') &
    (final_df['tier_liver'] != 'Tier_1_HighConfidence')
].head(15)
for _, row in top_t2.iterrows():
    report.append(f"- **{row['Symbol']}** ({row['functional_category'].replace('_',' ')}): master score {row['master_score']:.1f}")
report.append("")

report.append("## Suggested Next Steps")
report.append("1. **CRISPR/qPCR validation** of Immediate Priority genes (shared strong cis-eQTLs).")
report.append("2. **Allele-specific expression (ASE)** check from RNA-seq BAMs for top 10 genes.")
report.append("3. **Colocalization** with metabolic trait GWAS if available (e.g., COLOC for Fam111a, Scd3, Lcn2).")
report.append("4. **Gene set enrichment** on the full final_actionable_candidates.csv list using clusterProfiler/Enrichr.")
report.append("5. **Fine-mapping** with SuSiE on the top 5 loci to identify causal variants.")
report.append("6. **Trans-band follow-up**: For the top 3 cross-tissue trans-bands (Chr8:40-45 Mb, Chr5:160-165 Mb, Chr4:105-115 Mb), inspect raw read depth and founder haplotype effects to distinguish CNVs from mapping artifacts.")
report.append("7. **Pangenomic re-mapping**: Propose mapping RNA-seq reads from trans-band loci against a rat pangenome graph. If signal collapses, it is a mapping artifact; if it persists, it tags a genuine structural variant acting as a master regulator.")

with open('final_recommendations.md', 'w') as f:
    f.write('\n'.join(report))
print("Saved: final_recommendations.md")

print("\n=== Done ===")

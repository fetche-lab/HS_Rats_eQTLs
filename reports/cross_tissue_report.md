# Cross-Tissue eQTL Integration Report
Date: 2026-05-28

## 1. Dataset Overlap
- Adipose-only genes: 694
- Liver-only genes: 373
- Genes in both tissues: 197

## 2. Genomic Context Distribution

### Adipose
- proximal: 500 (56.1%)
- promoter: 193 (21.7%)
- trans: 101 (11.3%)
- distal_cis: 97 (10.9%)

### Liver
- proximal: 324 (56.8%)
- promoter: 128 (22.5%)
- distal_cis: 62 (10.9%)
- trans: 56 (9.8%)

## 3. Pleiotropy Classification (genes in both tissues)
- shared_cis_conserved: 135
- trans_both_hotspot: 33
- shared_cis_flip: 28
- tissue_specific_cis: 1

## 3b. Trans-Band Hotspots
Total trans-bands identified: 73
- Chr8:40-45 Mb: 3 genes, top LOC691532 (-logP=169.4) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE
- Chr5:160-165 Mb: 2 genes, top Eno1-ps20 (-logP=160.6) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; PSEUDOGENE_TRANS; LOC_UNCHARACTERIZED
- Chr15:50-55 Mb: 3 genes, top LOC120101677 (-logP=141.1) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE
- Chr7:130-135 Mb: 3 genes, top RGD1566212 (-logP=138.4) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE
- Chr1:210-215 Mb: 1 genes, top Rfk-ps3 (-logP=150.6) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; PSEUDOGENE_TRANS
- Chr12:5-10 Mb: 3 genes, top LOC108349532 (-logP=128.8) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; PSEUDOGENE_TRANS; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE
- Chr19:10-15 Mb: 4 genes, top LOC291863 (-logP=101.4) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE
- Chr7:0-5 Mb: 9 genes, top LOC680933 (-logP=50.0) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE; NEAR_TELOMERE
- Chr8:55-60 Mb: 1 genes, top LOC100361008 (-logP=126.9) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; LOC_UNCHARACTERIZED
- Chr5:130-135 Mb: 3 genes, top Gpbp1l2 (-logP=97.7) [CROSS-TISSUE]
  Flags: CROSS_TISSUE_TRANS; MEGA_EFFECT; PSEUDOGENE_TRANS; LOC_UNCHARACTERIZED; MASTER_REGULATOR_CANDIDATE

## 4. Shared Conserved cis-eQTLs (same peak, same direction)
Count: 135

| Symbol | Description | Adipose -logP | Adipose ES | Liver -logP | Liver ES | Peak Dist (Mb) |
|--------|-------------|---------------|------------|-------------|----------|----------------|
| Fam111a | FAM111 trypsin like peptidase A | 167.2 | 3.27 | 133.9 | 3.87 | 0.000 |
| RT1-N2 | RT1 class Ib, locus N2 | 118.2 | -2.91 | 133.4 | -3.82 | 0.000 |
| LOC103693260 | uncharacterized LOC103693260 | 153.2 | 3.56 | 129.8 | 3.28 | 0.000 |
| H2ac18 | histone cluster 2 H2A family member A2 | 92.6 | -2.85 | 98.2 | -3.49 | 0.000 |
| Ifit1 | interferon-induced protein with tetratricopeptide  | 150.0 | 2.34 | 129.2 | 3.20 | 0.000 |
| LOC102555965 | uncharacterized LOC102555965 | 115.0 | -1.87 | 142.9 | -4.92 | 0.000 |
| LOC120099939 | uncharacterized LOC120099939 | 119.5 | 2.05 | 175.6 | 4.16 | 0.000 |
| Lcn2 | lipocalin 2 | 120.7 | -2.19 | 89.9 | -2.44 | 0.000 |
| LOC100912642 | cytochrome P450 2J3-like | 134.1 | -2.02 | 59.2 | -1.17 | 0.000 |
| Mis18a | MIS18 kinetochore protein A | 123.8 | 1.83 | 113.8 | 2.36 | 0.000 |
| LOC120093570 | uncharacterized LOC120093570 | 144.8 | 2.73 | 31.3 | 0.76 | 0.000 |
| Hba-a3 | hemoglobin alpha, adult chain 3 | 85.8 | -1.95 | 99.0 | -2.48 | 0.000 |
| Ubd | ubiquitin like modifier D | 74.2 | 1.75 | 96.3 | 2.92 | 0.000 |
| LOC102550180 | uncharacterized LOC102550180 | 149.1 | 1.86 | 69.4 | 1.02 | 0.000 |
| Scd3 | stearoyl-coenzyme A desaturase 3 | 83.4 | 3.43 | 20.1 | 0.54 | 0.000 |
| Nat3 | N-acetyltransferase 3 | 93.3 | -1.74 | 70.1 | -1.79 | 0.000 |
| LOC120102367 | endogenous retrovirus group K member 21 Env polypr | 77.4 | -2.32 | 64.2 | -1.43 | 0.000 |
| Abca17 | ATP-binding cassette, subfamily A (ABC1), member 1 | 120.5 | -1.34 | 36.2 | -0.93 | 0.000 |
| RT1-CE12 | RT1 class I, locus CE12 | 137.2 | 1.14 | 119.9 | 2.60 | 0.000 |
| Galnt12 | polypeptide N-acetylgalactosaminyltransferase 12 | 109.8 | -1.84 | 52.5 | -1.08 | 0.000 |

## 5. Tissue-Specific Strong cis-eQTLs
Count: 1

| Symbol | Tissue | -logP | ES | Context | Description |
|--------|--------|-------|----|---------|-------------|
| LOC108351885 | Adipose | 44.7 | -0.76 | proximal | uncharacterized LOC108351885 |

## 6. Top Master Candidates (cross-tissue scored)
Total scored: 1231

| Rank | Symbol | Pleiotropy | Adipose Tier | Liver Tier | Master Score |
|------|--------|------------|--------------|------------|--------------|
| 1 | Fam111a | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 87.9 |
| 2 | RT1-N2 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 82.1 |
| 3 | LOC103693260 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 80.9 |
| 4 | H2ac18 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 79.1 |
| 5 | Ifit1 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 77.8 |
| 6 | LOC102555965 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 75.6 |
| 7 | LOC120099939 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 75.4 |
| 8 | Lcn2 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 73.7 |
| 9 | LOC100912642 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 71.1 |
| 10 | Mis18a | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 70.7 |
| 11 | LOC120093570 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 69.0 |
| 12 | Hba-a3 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 68.2 |
| 13 | Ubd | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 67.5 |
| 14 | LOC102550180 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 67.1 |
| 15 | Scd3 | shared_cis_conserved | Tier_1_HighConfidence | Tier_2_Validate | 67.1 |
| 16 | Nat3 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 67.0 |
| 17 | LOC120102367 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 66.4 |
| 18 | Abca17 | shared_cis_conserved | Tier_1_HighConfidence | Tier_2_Validate | 65.8 |
| 19 | RT1-CE12 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 65.7 |
| 20 | Galnt12 | shared_cis_conserved | Tier_1_HighConfidence | Tier_2_Validate | 65.7 |
| 21 | Syce1 | shared_cis_conserved | Tier_1_HighConfidence | Tier_1_HighConfidence | 65.5 |
| 22 | RGD1562844 | shared_cis_conserved | Tier_2_Validate | Tier_2_Validate | 65.5 |
| 23 | RGD1563402 | shared_cis_conserved | Tier_1_HighConfidence | Tier_2_Validate | 65.3 |
| 24 | Zfp94l1 | shared_cis_conserved | Tier_1_HighConfidence | Tier_2_Validate | 64.9 |
| 25 | Scrn1 | shared_cis_conserved | Tier_2_Validate | Tier_1_HighConfidence | 64.1 |
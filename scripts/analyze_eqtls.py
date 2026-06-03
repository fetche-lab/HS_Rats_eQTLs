import pandas as pd
import numpy as np

# Load files
adipose = pd.read_csv('HSNIH-Palmer_r4_HSNIH-Palmer_Adipose_RNA-Seq__Feb26__rlog_table.csv')
liver = pd.read_csv('HSNIH-Palmer_r4_HSNIH-Palmer_Liver_RNA-Seq__Feb26__rlog_table.csv')

for name, df in [('Adipose', adipose), ('Liver', liver)]:
    print(f"\n{'='*60}")
    print(f"Dataset: {name}")
    print(f"Shape: {df.shape}")
    print("\nColumns:", df.columns.tolist())
    print("\nFirst few rows:")
    print(df.head())
    print("\nNumeric summary:")
    print(df.describe())
    
    # Parse chromosomes and positions
    df['Gene_Chr'] = df['Location'].str.split(':').str[0]
    df['Gene_Pos'] = df['Location'].str.split(':').str[1].astype(float)
    df['Peak_Chr'] = df['Peak Location'].str.split(':').str[0]
    df['Peak_Pos'] = df['Peak Location'].str.split(':').str[1].astype(float)
    
    # Cis vs Trans: same chromosome and within 1 Mb or 4 Mb? Common definition for cis is within 1-4 Mb of gene body.
    # Let's use 4 Mb for a generous cis definition.
    df['Distance'] = np.where(df['Gene_Chr'] == df['Peak_Chr'], 
                              np.abs(df['Gene_Pos'] - df['Peak_Pos']), 
                              np.inf)
    df['cis_4mb'] = (df['Gene_Chr'] == df['Peak_Chr']) & (df['Distance'] <= 4.0)
    df['cis_1mb'] = (df['Gene_Chr'] == df['Peak_Chr']) & (df['Distance'] <= 1.0)
    
    print(f"\nCis/Trans counts (4 Mb threshold):")
    print(df['cis_4mb'].value_counts())
    print(f"\nCis/Trans counts (1 Mb threshold):")
    print(df['cis_1mb'].value_counts())
    
    # Highest -logP
    print(f"\nTop 10 by Peak -logP:")
    print(df.nlargest(10, 'Peak -logP')[['Symbol','Location','Peak -logP','Peak Location','Effect Size','cis_4mb','Distance']])
    
    # Correlation between -logP and Effect Size
    corr = df['Peak -logP'].corr(df['Effect Size'].abs())
    print(f"\nCorrelation between |Effect Size| and Peak -logP: {corr:.3f}")

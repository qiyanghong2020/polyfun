import subprocess
# 20231221
import pandas as pd

# Define the file paths
file_path = "/data/hongqy/workdir/FineMapping/COPD/after_dbsnp_anno/all.clump_results.20231221.clumps"
gwas_file = '/data/hongqy/workdir/FineMapping/COPD/COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.leftover.h19.header.1_22.txt'

# Read the parquet file using Pandas
df = pd.read_csv(gwas_file, sep='\t')
# 转换 'CHR' 列为数值型
df['#CHR'] = pd.to_numeric(df['#CHR'], errors='coerce')

# 转换 'BP' 列为数值型
df['POS'] = pd.to_numeric(df['POS'], errors='coerce')

mean_N_default = int(df['N'].mean())
# 删除 'CHR' 或 'BP' 列中包含 NaN 的行
df = df.dropna(subset=['#CHR', 'POS'])
#print(df)
# Read and process the clumps file
with open(file_path, 'r') as file:
    for line in file:
        # Split the line into columns
        columns = line.strip().split()
        if not columns[0].startswith("#"):
            CHROM, POS = columns[:2]
            CHROM = CHROM.strip()
            POS = POS.strip()

            # Ensure POS is an integer
            try:
                POS = int(POS)
            except ValueError:
                # If POS is not an integer, skip this line
                continue

            # Calculate Pstart and Pend
            Pstart = max(1, POS - 100000)
            Pend = POS + 100000

            # Filter dataframe for the specified region
            region_df = df[(df['#CHR'] == int(CHROM)) & (df['POS'].between(int(Pstart), int(Pend)))]

            # Save to text file
            output_filename = f'chr{CHROM}_{POS}.chr{CHROM}_{Pstart}_{Pend}.gwas.origin.txt'
            region_df.to_csv(output_filename, sep='\t', index=False)





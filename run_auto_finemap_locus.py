import subprocess
# 20231221, by hongqy
import os
import pandas as pd

# Define the file paths
file_path = "/data/hongqy/workdir/FineMapping/COPD/gwas_clump_finemap_all_update_rs/all.clump_results.20240110.clumps"
parquet_file = '/data/hongqy/workdir/FineMapping/COPD/gwas_clump_finemap_all_update_rs/COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.leftover.h19.header.txt.gwas2vcf.vcf.tsv.parquet'

# Read the parquet file using Pandas
df = pd.read_parquet(parquet_file)
# 转换 'CHR' 列为数值型
df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')

# 转换 'BP' 列为数值型
df['BP'] = pd.to_numeric(df['BP'], errors='coerce')

mean_N_default = int(df['N'].mean())
# 删除 'CHR' 或 'BP' 列中包含 NaN 的行
df = df.dropna(subset=['CHR', 'BP'])
#print(df)
# Read and process the clumps file
with open(file_path, 'r') as file:
    for line in file:
        # Split the line into columns
        columns = line.strip().split()
        if not columns[0].startswith("#"):
            CHROM, POS, ID = columns[:3]
            CHROM = CHROM.strip()
            POS = POS.strip()
            ID = ID.strip()

            # Calculate the mean of the 'N' column for the specified conditions
            mean_N = df.loc[(df['CHR'] == int(CHROM)) & (df['BP'] == int(POS)), 'N'].mean()

            # Check for missing values and set to default if necessary
            if pd.isna(mean_N):
                mean_N = mean_N_default
            else:
                mean_N = int(mean_N)  # or just keep it as a float if precision is important

            print('mean_N, CHROM, POS=', mean_N, CHROM, POS)
            # Ensure POS is an integer
            try:
                POS = int(POS)
            except ValueError:
                # If POS is not an integer, skip this line
                continue

            # Calculate P1 and P2 based on POS value
            if POS < 1000000:
                P1 = 1
                P2 = 3000001
            else:
                P1 = (POS // 1000000) * 1000000 + 1
                P2 = P1 + 3000000

            # Calculate Pstart and Pend
            Pstart = max(1, POS - 100000)
            Pend = POS + 100000

            # Filter dataframe for the specified region
            region_df = df[(df['CHR'] == int(CHROM)) & (df['BP'].between(int(Pstart), int(Pend)))]

            # Save to text file, easy for next step -> plot
            os.makedirs('output/', exist_ok=True)
            output_filename_prefix = f'chr{CHROM}_{POS}.{ID}.chr{CHROM}_{Pstart}_{Pend}.chr{CHROM}_{P1}_{P2}'
            region_df.to_csv(output_filename_prefix + '.gwas.txt', sep='\t', index=False)

            # Try to run different causal num
            causal_num = 5
            for causal_num in range(1, 6):
                if causal_num == 1:
                    ld = ""  # When max_num_causal=1 omit the flags --geno and --ld (we cannot use LD information )
                else:
                    ld = f"--ld /data/hongqy/database/polyfun/finemap/chr{CHROM}_{P1}_{P2}.dbsnp153 "
                # Prepare the command for each line
                command = f"""
                mkdir -p LD_cache
                mkdir -p output
    
                p='/data/hongqy/workdir/FineMapping/polyfun-master/'
                input='{parquet_file}'
                python $p/finemapper.py {ld} \\
                    --sumstats $input \\
                    --n {mean_N} \\
                    --chr {CHROM} \\
                    --start {Pstart} \\
                    --end {Pend} \\
                    --method finemap \\
                    --finemap-exe /data/hongqy/workdir/FineMapping/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 \\
                    --max-num-causal {causal_num} \\
                    --cache-dir LD_cache \\
                    --non-funct \\
                    --out output/{output_filename_prefix}.causal_{causal_num}.finemapper.gz \\
                    --allow-missing """

                with open(f'{output_filename_prefix}.causal_{causal_num}.finemapper.sh', 'w') as command_file:
                    command_file.write(command)
                # Run the command using subprocess
                process = subprocess.Popen(command, shell=True, executable='/bin/bash')
                process.wait()






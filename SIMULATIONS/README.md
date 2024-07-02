# Simulating MAGEPRO

loop_batch_1000genes.sh
- loops through various sample sizes and h2g values and calls sbatch batch_sim_1000genes.sh

batch_sim_1000genes.sh
- for each sample size and h2g setting, runs run_simulation_1000genes.sh

run_simulation_1000genes.sh
- simulated a random gene, and run the following four scripts to simulate MAGEPRO

1. Sim_geno.py
- simulate genotypes using LD matrix from 1kg

2. Sim_SumStats.py
- simulate eqtl summary statistic for EUR and AMR. These become our external summary statistics
- uses simulated gene expression with pre-set heritability and simualted genotypes

3. run_susie.R 
- run susie_rss() on simulated summary satistics 

4. Sim_Model.py
- simulate MAGEPRO = run susie on external datasets, ridge regression combination

## Many functions used in our simulations are from the twas_sim tool 
- https://github.com/mancusolab/twas_sim
> twas_sim, a Python-based tool for simulation and power analysis of transcriptome-wide association analysis. Xinran Wang, Zeyun Lu, Arjun Bhattacharya, Bogdan Pasaniuc, Nicholas Mancuso, twas_sim, a Python-based tool for simulation and power analysis of transcriptome-wide association analysis, Bioinformatics, 2023;

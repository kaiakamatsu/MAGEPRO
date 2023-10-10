# Simulating MAGEPRO

loop_batch_1000genes.sh
- loops through various sample sizes and h2g values

batch_sim_1000genes.sh
- batch run jobs

run_simulation_1000genes.sh
- simulated a random gene, create plink files from 1kg, run the following three python scripts to simulate MAGEPRO

Sim_geno.py
- simulate genotypes using LD matrix from 1kg

Sim_SumStats.py
- simulate eqtl summary statistic for one population 
- uses simulated gene expression of pre-set heritability and simualted genotypes

Sim_Model.py
- simulate MAGEPRO 
- uses simualted gene expression, simulated genotypes, simulated summary statistics

diff_LD 
- directory for similar files, but to simulate a situation where the causal snp is different in two populations 

## Many functions used in our simulations are from the twas_sim tool 
- https://github.com/mancusolab/twas_sim
> twas_sim, a Python-based tool for simulation and power analysis of transcriptome-wide association analysis. Xinran Wang, Zeyun Lu, Arjun Bhattacharya, Bogdan Pasaniuc, Nicholas Mancuso, twas_sim, a Python-based tool for simulation and power analysis of transcriptome-wide association analysis, Bioinformatics, 2023;

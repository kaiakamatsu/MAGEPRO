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

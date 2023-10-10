# loop through various sample sizes and heritability values and submit jobs 

afr_sizes=80,160,240,400,500
heritability=(0.03 0.05 0.1 0.2)

for h in "${heritability[@]}"; do
    	
    sbatch batch_sim_1000genes.sh $afr_sizes $h

done

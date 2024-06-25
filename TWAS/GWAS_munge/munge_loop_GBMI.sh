pop=$1

traits=(Asthma COPD Gout HF IPF Stroke VTE)

for c in "${traits[@]}"; do
	../../ldsc/munge_sumstats.py --sumstats ${pop}/filtered_hm3/${c}_${pop}_GBMI.txt --N-cas-col N_case --N-con-col N_ctrl --nstudy n_dataset --snp rsid --a1 REF --a2 ALT --p inv_var_meta_p --frq all_meta_AF --signed-sumstats inv_var_meta_beta,0 --out ${pop}/munged/${c}_${pop}_GBMI_munged
done

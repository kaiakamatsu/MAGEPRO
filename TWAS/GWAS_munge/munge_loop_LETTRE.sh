pop=$1

pheno=(BAS HGB MCHC MPV RBC EOS LYM MCV NEU RDW HCT MCH MON PLT WBC)

for c in "${pheno[@]}"; do
	../../ldsc/munge_sumstats.py \
  		--sumstats ../${pop}_BCX_LETTRELAB/filtered_raw/BCX2_${c}_${pop}_GWAMA.txt \
  		--N-col n_samples \
  		--snp rs_number \
  		--a1 reference_allele \
  		--a2 other_allele \
  		--p p-value \
  		--frq eaf \
  		--signed-sumstats z,0 \
  		--out ../${pop}_BCX_LETTRELAB/munged/BCX2_${c}_${pop}_GWAMA_munged
done

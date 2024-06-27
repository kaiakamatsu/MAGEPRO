pop=$1

pheno=(BAS HGB MCHC MPV RBC EOS LYM MCV NEU RDW HCT MCH MON PLT WBC)

for c in "${pheno[@]}"; do
	Rscript filter_sumstats_lexorder.R ../${pop}_BCX_LETTRELAB/raw/BCX2_${c}_${pop}_GWAMA.out.gz ../${pop}_BCX_LETTRELAB/filtered_raw/BCX2_${c}_${pop}_GWAMA.txt
done

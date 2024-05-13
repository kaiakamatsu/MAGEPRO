Rscript makePOS.R

phenos="BAS HGB MCHC MPV RBC EOS LYM MCV NEU RDW HCT MCH MON PLT WBC"

for pheno in $phenos; do
	echo $pheno
	sbatch runTWAS.sh $pheno
done

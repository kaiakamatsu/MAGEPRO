phenos="BAS HGB MCHC MPV RBC EOS LYM MCV NEU RDW HCT MCH MON PLT WBC"
#phenos="BAS"

for pheno in $phenos; do
	echo $pheno
	sbatch runTWAS.sh $pheno
done

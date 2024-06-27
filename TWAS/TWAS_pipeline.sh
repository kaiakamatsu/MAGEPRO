source ~/.bashrc
conda activate r_env

weights_dir=$1
genebed_dir=$2
letter_anc=$3
gbmi_anc=$4
weights_name=$5
#ld=${home}/GENOTYPES_GE_data/GENOA/GENOTYPES/phg000927.v1.NHLBI_GENOA_GWAS_AA.genotype-calls-matrixfmt.c1/AA_genotypes_combined/maf_hwe_rate_relatedness_HM3_filtered_people_with_GE_inGEU_YRI_redo/GENOA_AA_chr
ld=$6

Rscript list_wgt.R $weights_dir $weights_name

Rscript makePOS.R $genebed_dir $weights_name

phenos_letter="BAS HGB MCHC MPV RBC EOS LYM MCV NEU RDW HCT MCH MON PLT WBC"
phenos_gbmi="Asthma COPD Gout HF IPF Stroke VTE"

for pheno in $phenos_letter; do
	echo $pheno
	gwas=/expanse/lustre/projects/ddp412/kakamatsu/gwas/${letter_anc}_BCX_LETTRELAB/munged/BCX2_${pheno}_${letter_anc}_GWAMA_munged.sumstats.gz
	out_base=/expanse/lustre/projects/ddp412/kakamatsu/fusion_multipop/output_LETTER_${letter_anc}_${weights_name}
	sbatch runTWAS_SINGLE.sh $pheno $gwas $weights_dir ${out_base}_SINGLE $ld ${weights_name}
	sbatch runTWAS_SuSiE.sh $pheno $gwas $weights_dir ${out_base}_SuSiE $ld ${weights_name}
	sbatch runTWAS_PRSCSx.sh $pheno $gwas $weights_dir ${out_base}_PRSCSx $ld ${weights_name}
	sbatch runTWAS_MAGEPRO.sh $pheno $gwas $weights_dir ${out_base}_MAGEPRO $ld ${weights_name}
	sbatch runTWAS.sh $pheno $gwas $weights_dir ${out_base}_ALL $ld ${weights_name}
done


for pheno in $phenos_gbmi; do
	echo $pheno
	gwas=/expanse/lustre/projects/ddp412/kakamatsu/gwas/GBMI/${gbmi_anc}/munged/${pheno}_${gbmi_anc}_GBMI_munged.sumstats.gz
	out_base=/expanse/lustre/projects/ddp412/kakamatsu/fusion_multipop/output_GBMI_${gbmi_anc}_${weights_name}
	sbatch runTWAS_SINGLE.sh $pheno $gwas $weights_dir ${out_base}_SINGLE $ld ${weights_name}
	sbatch runTWAS_SuSiE.sh $pheno $gwas $weights_dir ${out_base}_SuSiE $ld ${weights_name}
	sbatch runTWAS_PRSCSx.sh $pheno $gwas $weights_dir ${out_base}_PRSCSx $ld ${weights_name}
	sbatch runTWAS_MAGEPRO.sh $pheno $gwas $weights_dir ${out_base}_MAGEPRO $ld ${weights_name}
	sbatch runTWAS.sh $pheno $gwas $weights_dir ${out_base}_ALL $ld ${weights_name}
done

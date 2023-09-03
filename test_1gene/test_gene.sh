gene=$1

row=$(cat /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/Genes_Expressed_in_Whole_Blood.txt | nl | grep $gene | awk '{print $1 + 1}')

zcat /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/Whole_Blood.v8.normalized_expression.bed.gz | head -n $row | tail -n 1 | cut -f 8,15,16,28,32,35,41,50,57,71,74,77,102,104,108,118,120,130,146,159,194,197,203,231,246,262,276,283,285,293,306,320,322,323,327,336,337,344,364,376,377,383,397,400,426,431,441,452,461,466,475,503,511,520,521,523,529,531,541,542,549,551,553,554,575,583,589,592,593,612,615,617,619,625,630,655,657,662,665,667 > gedonors.txt

/expanse/lustre/projects/ddp412/kakamatsu/plink2 --bfile /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/plink_cissnps_AFR/$gene --make-bed --keep /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/All_Individuals_Whole_Blood.txt --out /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene

rm /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene.log

paste --delimiters='\t' <(cut -f1-5 /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene.fam) <(cat gedonors.txt | sed 's/\t/\n/g') > /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene.mod.fam

mv /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene.mod.fam /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene.fam

Rscript ../ComputeMultiPopWeights.R --gene $gene --bfile /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/AFR_Whole_Blood/$gene --covar /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/Covar_All_Whole_Blood.txt --hsq_p 1 --tmp /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/tmp_Whole_Blood/"$gene"_Whole_Blood --out /expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/weights_trials/Whole_Blood/Whole_Blood.$gene --PATH_gcta /expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust --verbose 1 --PATH_plink ../../plink --sumstats_dir /expanse/lustre/projects/ddp412/kakamatsu/MAGEPRO_datasets --sumstats eqtlgen,genoa,mesahis,eurgtex,mesaafr,ota_CD16p_Mono,ota_CL_Mono,ota_LDG,ota_mDC,ota_Mem_CD4,ota_Mem_CD8,ota_Naive_B,ota_Naive_CD4,ota_Naive_CD8,ota_Neu,ota_NK,ota_pDC,ota_Plasmablast --models SINGLE,META,MAGEPRO --ss 31684,1031,352,574,233,416,416,416,416,416,416,416,416,416,416,416,416,416 --cell_type_meta ota 


#--crossval 1

#Create multipop weight files for GTEx v8 
#Begin an interactive session before running the script

#Set up r environment 
source ~/.bashrc
conda activate r_env

#Identify downsampled individuals and cis window 

Rscript 1_v8.NoDownsample.R #done (finds individuals with both genotype and gene expression data in GTEx)

#Only run once or everytime scratch directory deletes itself
run=2 #run: 1, don't run: 2 
if [ $run -eq 1 ]; 
then
 Rscript 2_v8.GetMeasuredGenes.R 
 #SNP extraction (tissue independent)
 ref_dir=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/data/GTEx_plink
 beddir=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/gene_beds/
 plinkdir=/expanse/lustre/scratch/kakamatsu/temp_project/GTExTEMP/plink_cissnps_AFR
 mkdir $plinkdir
 for i in $beddir/* #looping through every bed file created in Rscript 2
 do
 echo $i 
 name=`basename $i | sed 's/.bed//g'` #get gene name
 chrom=$(awk '{print $1}' $i) #get chromosome number
 plink_file=GTEx_v8_genotype_AFR_HM3_exclude_dups.${chrom} 
 ../../plink2 --bfile $ref_dir/$plink_file --extract bed1 $i --out $plinkdir/$name --make-bed #extract SNPs from GTEx around the region specified and create bed 
 #NOTE ABOUT PLINK --extract bed1 <file>: file has to be in this specific BED format 
 # chromosome code, start range, end range, set ID, group label 
 rm $plinkdir/$name.log
 done
fi

Rscript 3_v8.covar.AFR.R #done
Rscript assignBatch.R #done

multipopweights.sh

#after jobs run 
#Rscript processFlexPop.R 

#compute multipop weights in batches 
dataset=eqtlgen,ota,his,eur,mesa,genoa
batchfile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/genes_assign_Whole_Blood2.txt
batches=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $batchfile) #from 1-20 (batch number in the second column of the batchfile)
#for batch in $batches
#do
	#echo $batch
	#sbatch 2_job_benchmark.sh Whole_Blood $batch $dataset
#done


sbatch 2_job_benchmark.sh Whole_Blood 1 $dataset
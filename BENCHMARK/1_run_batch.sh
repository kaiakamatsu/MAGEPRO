#compute multipop weights in batches 
dataset=eqtlgen,ota,his,eur,mesa,genoa
batchfile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/gene_subset.txt
batches=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $batchfile) 
for batch in $batches
do
	echo $batch
	sbatch 2_job_benchmark.sh Whole_Blood $batch $dataset
done

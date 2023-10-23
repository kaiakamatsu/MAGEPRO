#compute multipop weights in batches 
dataset=eqtlgen,ota,his,eur,mesa,genoa
#batchfile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/genes_assign_Whole_Blood3.txt
batchfile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/benchmark/rerun_some/genes_assign_Whole_Blood.redo.txt
batches=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $batchfile) 
for batch in $batches
do
	if [ $batch -ne 1 ]; then
		echo $batch
		sbatch 2_job_benchmark.sh Whole_Blood $batch $dataset
	fi
done

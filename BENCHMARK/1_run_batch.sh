home=$1 #path to home dir
ldref=$2
gefile=$3
dataset=eqtlgen,his,mesa
batchfile=${home}/intermediate/Genes_Assigned.txt
batches=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $batchfile) 
for batch in $batches
do
	sbatch 2_job_benchmark.sh $batch $dataset $home $ldref $gefile
done

#sbatch 2_job_benchmark.sh 1 $dataset $home $ldref $gefile

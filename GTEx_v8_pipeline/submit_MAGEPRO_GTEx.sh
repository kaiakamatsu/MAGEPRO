#run jobs to compute MAGEPRO weights - in batches 

dataset=$1
batchfile=$2

ge=$3
scratch=$4 #make wd and temp, get plink
intermed=$5
output=$6
plink_exec=$7
plink_exec1=$8
gcta=$9

echo "RUNNING JOBS"
batches=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $batchfile)
for batch in $batches
do
	sbatch runMAGEPRO_GTEx.sh $batch $dataset $ge $scratch $intermed $output $plink_exec $plink_exec1 $gcta
done

#sbatch runmultipopweights.sh 1 $dataset $ge $scratch $intermed $output $plink_exec $plink_exec1 $gcta


#compute multipop weights in batches 

batchfile=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/intermedfiles/genes_assign_Whole_Blood2.txt
batches=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $batchfile) #from 1-20 (batch number in the second column of the batchfile)
for batch in $batches
do
	echo $batch
	sbatch runmultipopweights.sh Whole_Blood $batch
done


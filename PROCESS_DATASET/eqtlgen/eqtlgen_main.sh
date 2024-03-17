smr=$1
eqtl_smr=$2
out_chr=$3
snp=$4
out_gene=$5

#for chr in {1..22}
#do 
#	echo $chr
#	$smr/smr-1.3.1 --beqtl-summary $eqtl_smr/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --query 1 --probe-chr $chr --out $out_chr/$chr
#done

for i in {1..22}
do
	Rscript splitFilter.R $out_chr/$i.txt $snp $out_gene
done



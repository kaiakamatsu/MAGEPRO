# simulate a random gene and create plink files for that 1MB region 
numsafr_s=$1
IFS=',' read -ra numsafr <<< $numsafr_s
heritability=$2

randomgenes=random_genes_multi_h${heritability}
mkdir $randomgenes

#causal snp is the same 
for ((i=1; i<=1000; i++))
do 

	rm -rf ${randomgenes}/*

	# choose a random TSS and take +- 500 kb
	
	while true; do

		#random chr
		chr=$((1 + RANDOM % 22))
		echo $chr
	
		#pick a random variant and extract the position number
		file=/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/1000g/1000G_EUR/1000G_EUR_$chr.bim
		num_variants=$(wc -l < "${file}")
		line=$((1 + RANDOM % num_variants))
		random_variant=$(sed -n "${line}p" "${file}")
		random_position=$(echo "${random_variant}" | awk '{ print $4 }')

		start_pos=$(($random_position - 500000))
		if [ $start_pos -lt 0 ]; then 
			start_pos=0
		fi
		end_pos=$(($random_position + 500000))

		echo $start_pos
		echo $end_pos

		#EUR 
		../../../plink --bfile /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/1000g/1000G_EUR/1000G_EUR_$chr --chr $chr --from-bp $start_pos --to-bp $end_pos --make-bed --out ${randomgenes}/EUR_1KG_chr${chr}_${random_position} --extract hapmapsnps.txt

		#AFR
		../../../plink --bfile /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/1000g/1000G_AFR/1000G_AFR_$chr --chr $chr --from-bp $start_pos --to-bp $end_pos --make-bed --out ${randomgenes}/AFR_1KG_chr${chr}_${random_position} --extract hapmapsnps.txt

		if [ ! -e ${randomgenes}/EUR_1KG_chr${chr}_${random_position}.bim ] || [ ! -e ${randomgenes}/AFR_1KG_chr${chr}_${random_position}.bim ]; then
    			echo "no snps remaining, pick another gene"
			rm -rf ${randomgenes}/* 
			continue
		fi

		EUR_variants_random_gene=$(wc -l < "${randomgenes}/EUR_1KG_chr${chr}_${random_position}.bim")
		AFR_variants_random_gene=$(wc -l < "${randomgenes}/AFR_1KG_chr${chr}_${random_position}.bim")

		if [ $EUR_variants_random_gene -lt 50 ]; then
			echo "not enough snps in chosen gene, pick another gene"
			rm -rf ${randomgenes}/*
			continue
		elif [ $AFR_variants_random_gene -lt 50 ]; then
			echo "not enough snps in chosen gene, pick another gene"
			rm -rf ${randomgenes}/*
			continue
		else
			echo "successfully picked a random gene"
			break
		fi
	done	
	
	min_snps=$((EUR_variants_random_gene < AFR_variants_random_gene ? EUR_variants_random_gene : AFR_variants_random_gene))
	#echo $min_snps
	while true; do
    		index_causal1=$((RANDOM % $min_snps))
    		index_causal2=$((RANDOM % $min_snps))
    		if [ "$index_causal1" -ne "$index_causal2" ]; then
        		break
    		fi
	done
	causal_indices="$index_causal1,$index_causal2"
	echo $causal_indices

	for numafr in "${numsafr[@]}"; do

		simgeno=simulated_genotypes_multi_${numafr}_h${heritability}

		sumstatsdir=sumstats_multi_${numafr}_h${heritability}

		tempdir=temp_multi_${numafr}_h${heritability}

		if [ $i -eq 1 ]; then
			mkdir $simgeno
			mkdir $sumstatsdir
			mkdir $tempdir
		else 
			rm -rf ${simgeno}/*
			rm -rf ${sumstatsdir}/*
			rm -rf ${tempdir}/*
		fi

		# simulate genotypes using python script
		python3 Sim_geno.py ${randomgenes}/EUR_1KG_chr${chr}_${random_position} 500 EUR ${simgeno}
		python3 Sim_geno.py ${randomgenes}/AFR_1KG_chr${chr}_${random_position} $numafr AFR ${simgeno}

		# for both AFR and EUR with the same causal variant for now

		# simulate eqtl analysis to produce EUR full sum stats 
		#python3 Sim_SumStats.py <path to plink files> <sample size> <population> <simulated genotypes> <number of causal snps> <index of causal> <simulation number> <heritability of gene> <sample size of target pop> <sumstats dir> <tempdir>
		python3 Sim_SumStats_multi.py ${randomgenes}/EUR_1KG_chr${chr}_${random_position} 500 EUR ${simgeno}/simulated_genotypes_EUR.csv 2 $causal_indices $i $heritability $numafr $sumstatsdir $tempdir

		# simulate AFR gene models with lasso 
		# use MAGEPRO 
		#python3 Sim_Model.py <path to plink files> <sample size> <population> <simulated genotypes> <sumstats> <number of causal snps> <index of causal> <simulation number> <heritability of gene>
		python3 Sim_Model_multi.py ${randomgenes}/AFR_1KG_chr${chr}_${random_position} $numafr AFR ${simgeno}/simulated_genotypes_AFR.csv ${sumstatsdir}/sumstats_EUR.csv 2 $causal_indices $i $heritability $tempdir
	done
done

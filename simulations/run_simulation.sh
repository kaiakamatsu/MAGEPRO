# simulate a random gene and create plink files for that 1MB region 
run=0

rm -rf random_genes/*

if [ $run -eq 0 ]; then
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
		../../../../plink --bfile /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/1000g/1000G_EUR/1000G_EUR_$chr --chr $chr --from-bp $start_pos --to-bp $end_pos --make-bed --out random_genes/EUR_1KG_chr${chr}_${random_position}

		#AFR
		../../../../plink --bfile /expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/1000g/1000G_AFR/1000G_AFR_$chr --chr $chr --from-bp $start_pos --to-bp $end_pos --make-bed --out random_genes/AFR_1KG_chr${chr}_${random_position}
	
		EUR_variants_random_gene=$(wc -l < "random_genes/EUR_1KG_chr${chr}_${random_position}.bim")
		AFR_variants_random_gene=$(wc -l < "random_genes/AFR_1KG_chr${chr}_${random_position}.bim")

		if [ $EUR_variants_random_gene -lt 50 ]; then
			rm -rf random/genes/*
			continue
		elif [ $AFR_variants_random_gene -lt 50 ]; then
			rm -rf random/genes/*
			continue
		else
			echo "successfully picked a random gene"
			break
		fi
	done
	
fi


# simulate genotypes using python script 
python3 Sim_geno.py 'random_genes/EUR_1KG_chr${chr}_${random_position}' 1500 'EUR'
python3 Sim_geno.py 'random_genes/AFR_1KG_chr${chr}_${random_position}' 1500 'AFR'

# simulate total expression (with h2 = 0.05 betas and simulated genotypes)

	# for both AFR and EUR with the same causal variant

	# for both AFR and EUR with different causal variant


# simulate eqtl analysis to produce EUR full sum stats 



# simulate AFR gene models with lasso 



# use MAGEPRO 





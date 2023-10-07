# Extracting and preparing data from eqtlgen

Download eqtlgen data from ...
- https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz

Acquire smr software tool executable from ...
- https://yanglab.westlake.edu.cn/software/smr/#Download

We have prepared a script to extract the eqtl summary statistics and create gene-specific files

bash eqtlgen_main.sh 
  - path to smr software executable
  - path to directory containing eqtl sumstat file
  - path to output intermediate files: data split by chr 
  - path to snp reference file: SNPs to keep in analysis 
  - path to output gene specific files


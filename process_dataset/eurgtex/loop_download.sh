# download and configure gsutil from https://cloud.google.com/storage/docs/gsutil_install#linux

#gsutil -u <google cloud project> cp gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Whole_Blood.v8.EUR.allpairs.chr$c.parquet <output dir>
output=$1

for c in {1..22}
do 
	echo $c
	gsutil -u labproject-394420 cp gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Whole_Blood.v8.EUR.allpairs.chr$c.parquet $output # REPLACE PROJECT WITH YOUR GOOGLE CLOUD PROJECT	
done

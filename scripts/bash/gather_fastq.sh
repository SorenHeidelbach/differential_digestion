


fastq_dir="/user_data/sh/differential_digestion/data/raw/2023-02-08_AD-differential-digest/AD1_MboI/20230208_1510_X4_FAV45808_cb871377/fastq_pass"
fastq_out="/user_data/sh/differential_digestion/data/processed/NP/NP_MboI_R1041.fastq"

cat $fastq_dir/* > $fastq_out
if echo $(file $fastq_out) | grep -q "gzip"
then 
    mv $fastq_out $fastq_out.gz 
    gunzip $fastq_out.gz
fi

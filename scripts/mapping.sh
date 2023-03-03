WD="/user_data/sh/differential_digestion"
cd $WD

dir_fastq="data/processed/NP"
dir_out="mmlong2/results/mapping"
assembly="/user_data/sh/differential_digestion/mmlong2/tmp/polishing/asm_pol_lenfilt.fasta"
reads_to_map=$(ls $dir_fastq | grep MboI)

module load minimap2/2.24-GCCcore-10.2.0
mkdir -p $dir_out
for reads in $reads_to_map; do
    sam="$dir_out/${reads%.*}.sam"
    fastq=$dir_fastq/$reads
    echo "mapping $fastq to $assembly and outputting to $sam"
    minimap2 -ax map-ont -t 100 $assembly $fastq > $sam
done
module purge


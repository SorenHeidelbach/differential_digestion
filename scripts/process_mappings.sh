



wd=/user_data/sh/differential_digestion/mmlong2/results/mapping
mappings=$(ls $wd | grep MboI)
module load SAMtools/1.14-GCC-10.2.0
for mapping in $mappings; do
    echo $mapping
    mapping=$(basename $mapping)
    mapping=${mapping%.sam}
    #echo "sorting"
    #samtools sort -t 50 $wd/$mapping.sam > $wd/$mapping.sorted.sam
    echo "converting to bam"
    samtools view -t 50 -bS $wd/$mapping.sorted.sam > $wd/$mapping.sorted.bam
    echo "indexing"
    samtools index $wd/$mapping.sorted.bam
    echo "coverage"
    samtools coverage $wd/$mapping.sorted.bam -o $wd/$mapping.cov.tsv
done
module purge

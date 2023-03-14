import os
import re
from pathlib import Path
configfile: "config.yaml"

threads = config["threads"]
assm = config["assembly"]

np_reads = os.listdir(config["NP_reads"])
il_reads = os.listdir(config["IL_reads"])
pb_reads = os.listdir(config["PB_reads"])
coverage_reads = il_reads + pb_reads + np_reads

reads_to_map = [Path(sample).stem for sample in coverage_reads if re.search('MboI|260|400bps_sub\.', sample)]
print(reads_to_map)

rule all:
    input:
        expand("data/processed/mappings/{sample}/depth.txt", sample=reads_to_map)

rule read_mapping:
    input:
        assm=assm,
        reads=config["NP_reads"] + "/{sample}.fastq"
    output: "data/processed/mappings/{sample}/align.sorted.bam"
    threads: threads
    shell:
        """
        module load minimap2/2.24-GCCcore-10.2.0
        module load SAMtools/1.14-GCC-10.2.0

        minimap2 -ax map-ont -t {threads} {input.assm} {input.reads} |
        samtools sort -t {threads} |
        samtools view -t {threads} -bS  > {output}
        samtools index {output}

        module purge
        """

rule mapping_depth:
    input:
        "data/processed/mappings/{sample}/align.sorted.bam"
    output:
        "data/processed/mappings/{sample}/depth.txt"
    threads: threads
    shell:
        """
        module load SAMtools/1.14-GCC-10.2.0

        samtools depth --threads {threads} -aa {input} > {output}

        module purge
        """
import os
import re
configfile: "config.yaml"

def listdir_fullpath(dir):
    return [os.path.join(dir, i) for i in os.listdir(dir)]

np_reads = listdir_fullpath(config["NP_reads"])
il_reads = listdir_fullpath(config["IL_reads"])
pb_reads = listdir_fullpath(config["PB_reads"])
polish_reads = [reads for reads in np_reads if "R1041_400bps" in reads]
coverage_reads = il_reads + pb_reads + np_reads


rule all:
    input:
        "results/coverage"

rule assembly:
    input:
        np_reads
    output:
        directory("results/Flye")
    shell:
        """
        module load Flye
        flye --nano-raw {input} --meta --out-dir {output} --threads 50
        module purge
        """

rule medaka:
    input:
        assm="results/Flye/assembly.fasta",
        reads=polish_reads
    output:
        "results/medaka"
    conda:
        "/user_data/sh/conda_env/medaka"
    shell:
        """
        medaka_consensus -i {input.reads} -d {input.assm} -o {output} -t 50 -m r1041_e82_400bps_sup_g615
        """

rule assembly_QC:
    input:
        rules.medaka.output
    output: 
        directory("results/QUAST")
    shell:
        """
        module load QUAST/4.6.3-foss-2018a-Python-3.6.4
        quast.py -o {output} {input}/consensus.fasta
        module purge
        """

rule read_mapping_np:
    input:
        assm=rules.medaka.output,
        reads=np_reads
    output:
        directory("results/read_mappings_np")
    shell:
        """
        module load minimap2/2.24-GCCcore-10.2.0
        for read in {input.reads}
        do
            name=$(basename -s .fastq $read)
            minimap2 -ax map-ont -t 50 {input.assm}/consensus.fasta $read | samtools sort -o {output}/$name.bam
        done
        module purge
        """
rule read_mapping_pb:
    input:
        assm=rules.medaka.output,
        reads=pb_reads
    output:
        directory("results/read_mappings_pb")
    shell:
        """
        module load minimap2/2.24-GCCcore-10.2.0
        for read in {input.reads}
        do
            name=$(basename -s .fastq $read)
            minimap2 -ax map-hifi -t 50 {input.assm}/consensus.fasta $read | samtools sort -o {output}/$name.bam
        done
        module purge
        """
rule read_mapping_il:
    input:
        assm=rules.medaka.output,
        reads=il_reads
    output:
        directory("results/read_mappings_il")
    shell:
        """
        module load Bowtie2/2.4.2-foss-2020b
        for read in {input.reads}
        do
            name=$(basename -s .fastq $read)
            bowtie2-build --threads 50 {input.assm}/consensus.fasta {input.assm}/index
            bowtie2 -p 50 -x {input.assm}/index -U $read | samtools sort -o {output}/$name.bam
        done
        module purge
        """

rule calculate_coverage:
    input:
        "results/read_mappings"
    output:
        directory("results/coverage")
    shell:
        """
        module load CoverM/0.6.0-foss-2020b
        coverm contig -b {input}/* -o {output}
        module purge
        """
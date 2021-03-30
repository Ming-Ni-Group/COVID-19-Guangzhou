#!/bin/bash
###################################
# This script is a snake pipeline that is used for  primer sequences trimming in the sequencing sequence based on the bed files corresponding to the primers used in the nCov amplification and genome assembly  based on a reference genome.
################################################################################################################################
#
# Examples:
#	 
#	    snakemake -s nCov_assembly_pipeline.py  -p   
#		
#	
################################################################################################################################

#For different samples, you need to change the following "REP_INDEX" ,it is a list of all sample names.
#For example:
REP_INDEX = ['ID001','ID002']


# You need to change the output path for each processing step.
rule all:
	input:
	## bwa to ref  genome
		expand("process/bwa/{samp}/{samp}.sorted.bam",samp=REP_INDEX),
	##bam_index
		expand("process/bwa/{samp}/{samp}.sorted.bam.bai",samp=REP_INDEX),
	##refGonm_faidx
		expand("ref/MN908947.3_genome.fna/MN908947.3_genome.fna.fai"),
	##ivarTrim_Primer_sort
		expand("process/ivar_trim_Primer/{samp}.primertrim.sorted.bam",samp=REP_INDEX),
	##ivarTrim_Primer_sortindex
		expand("process/ivar_trim_Primer/{samp}.primertrim.sorted.bam.bai",samp=REP_INDEX),
	##ivarTrim_Primer_sort_depth
		expand("process/ivartrimPrimer_depth/{samp}.primertrim.sorted.depth.txt",samp=REP_INDEX),
	##mpileup_consensus
		expand("process/mpileup_consensus/{samp}.consensus.fa",samp=REP_INDEX),	
	##ivarTrim_varscan2Indel
		expand("process/varscan2/{samp}/{samp}.ivar_trim.varscan2.mpileup2Indel.vcf",samp=REP_INDEX),

rule bwa_index:
	input:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna"
	output:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna.bwt"
	shell:
		"bwa  index  {input}"

rule bwa_run:
	input:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna",
		"data/{samp}/{samp}_350.fq1.gz",
		"data/{samp}/{samp}_350.fq2.gz",
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna.bwt"
	output:
		"process/bwa/{samp}/{samp}.sorted.bam"
	log:
		"process/bwa/{samp}/{samp}.log"
	shell:
		"bwa  mem  -t  8  {input[0]}  {input[1]}  {input[2]} \
		| samtools sort -@  8 |  samtools view -F 4 -o  {output}   2>{log}  "
 
rule bam_index:
	input:
		"process/bwa/{samp}/{samp}.sorted.bam"
	output:
		"process/bwa/{samp}/{samp}.sorted.bam.bai"
	shell:
		"samtools index  {input}"

rule refGonm_faidx:
	input:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna"
	output:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna.fai"
	shell:
		"samtools  faidx  {input}"

rule ivar_trim_Primer:
	input:
		"process/bwa/{samp}/{samp}.sorted.bam",
		"yinwubed/ncov_1k_kuo.bed",
	output:
		temp("process/ivar_trim_Primer/{samp}.primertrim.bam"),
		"process/ivar_trim_Primer/log/{samp}.primertrim.bam.log",
	params:
		"process/ivar_trim_Primer/{samp}.primertrim",
	shell:
		"ivar trim -e -i {input[0]} -b {input[1]} -p {params[0]}  1>{output[1]}"

	#ivar trim
	#Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]
	#Input Options    Description
           	#-i    (Required) Sorted bam file, with aligned reads, to trim primers and quality
           	#-b    (Required) BED file with primer sequences and positions
           	#-m    Minimum length of read to retain after trimming (Default: 30)
           	#-q    Minimum quality threshold for sliding window to pass (Default: 20)
           	#-s    Width of sliding window (Default: 4)
           	#-e    Include reads with no primers. By default, reads with no primers are excluded
	#Output Options   Description
           	#-p    (Required) Prefix for the output BAM file


rule ivarTrim_Primer_sort:
	input:
		"process/ivar_trim_Primer/{samp}.primertrim.bam",
	output:
		"process/ivar_trim_Primer/{samp}.primertrim.sorted.bam"
	shell:
		"samtools sort -@ 6 {input} -o {output}"

rule ivarTrim_Primer_sortindex:
	input:
		"process/ivar_trim_Primer/{samp}.primertrim.sorted.bam",
	output:
		"process/ivar_trim_Primer/{samp}.primertrim.sorted.bam.bai"
	shell:
		"samtools index  {input}"


rule ivarTrim_Primer_depth:
	input:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna",
		"process/ivar_trim_Primer/{samp}.primertrim.sorted.bam",
	output:
		"process/ivartrimPrimer_depth/{samp}.primertrim.sorted.depth.txt"
	shell:
		"samtools depth -a  -d 20000  --reference {input[0]} {input[1]}  >{output}"

 

rule mpileup_consensus:
	input:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna",
		"process/ivar_trim_Primer/{samp}.primertrim.sorted.bam",
	output:
		"process/mpileup_consensus/{samp}.consensus.fa"
	params:
		"process/mpileup_consensus/{samp}.consensus",
	shell:
		"samtools mpileup -A -d 20000 -B -Q 0 --reference  {input[0]}  {input[1]} | ivar consensus -p {params} -n N"

rule mpileup:
	input:
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna",
		"process/ivar_trim_Primer/sorted/{samp}.primertrim.sorted.bam",
		"ref/MN908947.3_genome.fna/MN908947.3_genome.fna.fai",
		"process/ivar_trim_Primer/sorted/{samp}.primertrim.sorted.bam.bai"
	output:
		"process/ivar_trim_Mpileup/{samp}.ivar_trim.mpileup"
	log:
		"process/ivar_trim_Mpileup/{samp}.ivar_trim.mpileup.log"
	shell:
		"samtools mpileup -A  -d 10000  -B -Q 0 --reference  {input[0]}  {input[1]}  1>{output} 2>{log}"     


rule varscan2_mpileup2Indel:
	input:
		"process/ivar_trim_Mpileup/{samp}.ivar_trim.mpileup"
	output:
		"process/varscan2/{samp}/{samp}.ivar_trim.varscan2.mpileup2Indel.vcf"
	shell:
		"varscan   pileup2indel  {input}  --min-coverage 50  --min-reads2  5    --min-var-freq  0.3   --variants  indel --output-vcf 1  >{output} "     


# Genomic-elucidation-of-a-COVID-19-resurgence-and-local-transmission-of-SARS-CoV-2-in-Guangzhou
## Summary
While China experienced a dramatic decline of coronavirus disease 2019 (COVID-19) at the outset of 2020, regional outbreaks continuously emerged in recent months. The resurgence of COVID-19 has also been seen in many other countries. In Guangzhou, a small outbreak emerged in March and April involving less than 100 residents, and a comprehensive and near-real-time genomic surveillance of SARS-CoV-2 was conducted. When confirmed cases among overseas travelers increased, public health measures were enhanced as shifting self-quarantine to central quarantine and SARS-CoV-2 testing for all overseas travelers. From 109 imported cases we found diverse viral variants distributing in the global viral phylogeny, which were usually shared within households but not among passengers on the same flight. The viral variants from travelers also provided the information of COVID-19 in regions where viral genomic surveillance were lacked. In contrast to the viral diversity of imported cases, local transmission was predominately attributed to two specific variants imported from Africa, including the local cases who reported no direct/indirect contact with imported cases. The introducing events of the virus were identified or deduced before enhanced measures were taken. These results show that the interventions were effective in containing the spread of SARS-CoV-2, and also ruled out the possibility of cryptic transmission of viral variants from the first wave in January and February. Our study provides evidence and emphasizes the importance of controls for oversea travelers in the context of the pandemic, and exemplifies how viral genomic data facilitates COVID-19 surveillance and prevention.

## Dependencies
MAFFT 7.458
IQ-Tree2 rc2
RAxML 8.2.12
ggtree 1.14.6
pangolin V1.1.8(https://github.com/hCoV-2019/pangolin)
Biopython
python 2.7
python 3.6
R version 3.5.1

## Repo Contents
data: preprocessed data used for analysis.
result: result files for figure visualization.
scripts: python script and R script for analysis and visulization.

## Phylogenetic analysis
1. Multi-Alignment of the SARS-CoV-2 genomes : `mafft --thread 10 --auto $fasta_file > $align_result_file`

2. Trim Alignment file: `python ./scripts/Trim_UTR_multialign_fasta.py $align_result_file $align_trim_file`

3. Iqtree phylogenetic analysis: `iqtree -s $align_trim_file --prefix iqtree_result -m GTR+F+R2 -B 1000 -cmax 30 -redo -T 24`

4. RAXML phylogenetic analysis: `raxmlHPC-PTHREADS -s $align_trim_file -n raxmltree_result -m GTRGAMMAI -f a -x 12345 -N 1000 -p 123456 -T 24 -k`

5. pangolin lineage analysis: `pangolin $fasta_file -t 24`


## Genome assembly and variations calling

1. To obtain consensus sequences and deletion mutations: `snakemake -s GenomeAssembly_And_indelCalling.pipeline.py  -p`

2. Calling iSNV and SNP mutations by using script "iSNV_calling.sh" according to https://github.com/generality/iSNV-calling: `bash iSNV_calling.sh  3,4,5`

## variations annotation

1. To convert iSNV table into vcf format(VCFv4.1) that can be recognized  by SnpEff software: `python ./scripts/Variantions_annotion/iSNVTable_2_vcf_allsample.py  -i  ./data/   -o $allsamplesVcfPath   -r   MN908947.3`
2. Annoting the mutation by using SnpEff software ,for example: `java -jar snpEff.jar ann  MN908947   $$allsamplesVcfPath/${sample}.vcf    >$outSnpeff/{sample}.snpeffAnno.vcf`
3. To Integrate all the samples' mutations annotion files into a single file: `python Snpeff_results_integrateTo1file -i   $outSnpeff/   -o  ./data/allsampMutatAnno.txt`

4. To annotate the indel that were filted:`java -jar snpEff.jar ann  MN908947   ${indelFPath}/${sample}.indel.vcf    >.${indelAnnoPath}/{sample}.indel.snpeffAnno.vcf`


## Citation
If you use data, results or conclusion from this work, please cite:

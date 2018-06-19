# Investigation of Developmental and Shh Gene Regulatory Regions: A Differential Expression Analysis of RNA+ChIP-Sequencing Data
##### Mauro Chavez
##### 5/1/2018

## Abstract

This lab involved the analysis of RNA-Seq and ChIP-Seq data from three different tissue types from a developing mouse. Differential expression analysis was used to identify genes that behaved differently across tissue types. This information was considered along side locations of histone modifications to identify their effect in gene regulation. The Shh gene and its ZRS enhancer region were highlighted as regions of interest. Upon analysis of ZRS region in different species, those without limbs contained altered ZRS sequence. Differential expression data is useful not only in this context, but also has machine learning applications to cancer research.

## Introduction

This week’s lab focused on studying the regulation of genes in early mouse development to understand how an embryo differentiates cell types during growth. More specifically, which genes and regulatory regions control the development of limbs like legs and arms vs. other tissues like those found in the brain. Once these regions were identified, they would be used to study snake genomes to uncover sequence changes potentially resulting in a loss of limbs.

Analysis was completed using RNA-Sequencing data that was captured by reverse transcription of mRNA into complementary cDNA. cDNA is sequenced and mapped to a reference transcriptome. Alternatively, a given transcriptome could have been split into kmers where kmer counting in reads could quantify gene expression levels. Genes with more reads mapped to them are understood to have higher expression than those with fewer mapped reads. From here statistics can be gathered on number of reads, fragments, or transcripts for each gene per million reads. Given this information for various genes and tissue types, differential expression analysis would identify genes with different behavior in different tissues. Highly conserved differentially expressed genes in limb and non-limb tissues were hoped to provide hints at the genes causing loss of limbs in snakes. To add extra data, ChIP sequencing information of H3K27Ac and H3K4me1 histones were included in the analysis to determine their contributions to gene regulation and appearance in different tissue types. 

This type of work is important because differential expression analysis of RNA-Sequencing data can be used in many different contexts. For example, in the paper “RNA-Seq Accurately Identifies Cancer Biomarker Signatures to Distinguish Tissue of Origin”, Wei et al. used RNA-Seq data to generate expression profiles for various cancer types. A machine learning approach was then used to build a classifier to identify tissue of origin for metastatic cancer cells given RNA-Seq data with 90.5% accuracy [7]. This strengthens the usefulness of RNA-Seq analysis as a tool that can be used to uncover the defining characteristics between groups of study. 


## Methods 

### RNA-seq datasets
Samples from hind limb (HL), fore limb (FL) and mid brain(MB) were taken from a developing mouse at day E10.5. For each tissue, two replicates were collected, sequenced, and labeled `Rep1` and `Rep2`. RNA-Seq paired end reads for each tissue type were split into separate files. Total RNA-sequencing data from these tissue types were down sampled to only contain data from chromosome 5. Each tissue type had read lengths of 50 bp with various numbers of reads ranging from 2.9 to 3.9 million paired end reads. Read length statistics are reported in Table 1. FASTQC [1] was run on all foreword and reverse reads for each replicate from each tissue type. In each case, all modules passed except Per base sequence content. There were nucleotide biases such that the difference between any two nucleotide’s frequencies exceeded 20% at a certain position. In each read file, biases leveled out at the 16th base pair. This is normal for RNA-Seq libraries due to overrepresented sequences or biases in random hexamer priming. GRCm38 also known as mm10 mouse reference genome was used in all analysis steps. 

| experiment  | read length | number of reads | number of paired end reads |
| - | ----------- | --------------- | -------------------------- |
| FL_Rep1_chr5 | 50 | 6358798 | 3179399 |
| FL_Rep2_chr5 | 50 | 5930058 | 2965029 |
| HL_Rep1_chr5 | 50 | 7865676 | 3932838 |
| HL_Rep2_chr5 | 50 | 5623826 | 2811913 |
| MB_Rep1_chr5 | 50 | 6509950 | 3254975 |
| MB_Rep2_chr5 | 50 | 6827878 | 3413939 | 

*Table 1: read length, number of reads, number of paired end reads for each tissue type and replicant*

### RNA-seq analysis 
STAR was run using 4 threads to align RNA-Seq reads to reference mm10 transcriptome. Kallisto[6] version 0.44.0 was then used to qualify transcript abundance quickly, turning fastq files into transcripts per million reads values. Script run_kallisto.sh ran the quantification algorithm using three threads for each tissue replicate. Kallisto quant was run using 100 bootstrap samples instead of the default 0 in order to increase the confidence in reported accuracy of gene expression analysis. The --gtf flag was used so that transcripts could be translated into genomic coordinates.  Reference genome mm10 in kallisto index format was provided to the –i flag for pseudoalignment. Both foreword and reverse paired end reads for each tissue and replicate was provided to kallisto. The output of our kallisto analysis was used in sleuth[2] to identify differentially expressed genes. Sleuth was run in an R environment, with paths set to kallisto results, metadata loaded, sleuth object created with `extra_bootstrap_summary` parameter set to true, and false discovery rate threshold set to less than or equal to 0.05. Significant results were outputted to a .tab file. 


### Enhancer analyses
The Integrative Genomics Viewer (IGV)[3] was used to visualize .tdf files with count per position information that were created by aligning RNA-Seq reads to reference transcriptome mm10 using STAR[8] with results processed using igvtools. IGV was launched with 750MB setting using the mm10 reference genome. .tdf files for each tissue and replicate were loaded into IGV in addition to H3K27ac and H3K4me1 ChIP-sequencing data for each tissue type . A PhyloP track for the mm10 chromosome 5 was also loaded in to compare sequence conservation across species. Multiple sequence alignment of ZRS enhancer sequences was accomplished using mafft[4]. mafft was run with --auto flag on .fa file of ZRS enhancer sequences. From there, web version 1.54 of Mview[5] was used to visualize alignment. Alignments were pasted into text box with DNA option selected. From there input parameters were set to automatic and output parameters were set to default. 

## Results

### RNA-seq analysis

pairwisePearson.sh was used to calculate Pearson correlation between TPM values across all pairwise kallisto results. The data from these calculations is presented in table 2.  RNA-Seq replicates of the same tissue type were always extremely comparable with correlation of 0.99+. Fore limb and hind limb tissues were most similar when comparing correlation of different tissue types. However, fore limb and mid brain were only 1% less correlated than fore limb was to hind limb in some cases. Mid brain and fore limb tissues were most dissimilar.

|       | FL_Rep1 | FL_Rep2 | HL_Rep1 | HL_Rep2 | MB_Rep1 | MB_Rep 2 |
| ----- | ------- | ------- | ------- | ------- | ------- | -------- |
| FL_Rep1 | 1|0.9971084853118|0.96069401846103|0.95701275529406|0.957516800544|0.94686987513034 |
| FL_Rep2 |0.9971084853118|1|0.95954358025087|0.95686850591509|0.96011384025694|0.95201712558146|
| HL_Rep1 |0.96069401846103|0.95954358025087|1|0.99561458942698|0.93632509599264|0.92713427649162|
| HL_Rep2 |0.95701275529406|0.95686850591509|0.99561458942698|1|0.93527689258929|0.92697081070712|
| MB_Rep1 |0.957516800544|0.96011384025694|0.93632509599264|0.93527689258929|1|0.99164218909938|
| MB_Rep2 |0.94686987513034|0.95201712558146|0.92713427649162|0.92697081070712|0.99164218909938|1|

*Table 2: Pairwise Pearson correlation of kallisto results*

After identifying differentially expressed genes with the help of sleuth, 1943 significant transcripts were identified.  Using IGV, expression across tissues and replicates were analyzed. In the case of Shh (see figure 1), it was highly expressed in the mid brain, and less expressed in hind limb and fore limb. The same story was also true for the Uchl1 gene. Expression levels for Uchl1 are reported in figure 2. Although in this case, the difference in level of expression in fore limb vs hind limb was slightly larger. This was also the case for Sparcl1 expression as shown in figure 3.  The Ubl3 gene was expressed very similarly across all tissue types as shown in figure 4. Expression for Nom1 shown in figure 5 is one gene that had higher expression in hind limb and fore limb. 

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/ShhExpression.png)

*Figure 1: Expression of Shh gene across all tissue types and replicates*

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/Uchl1Expression.png)

*Figure 2: Expression of Uchl1 gene across all tissue types and replicates*

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/Sparcl1Expression.png)

*Figure 3: Expression of Sparcl1 gene across all tissue types and replicates*

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/Ubl3Expression.png)

*Figure 4: Expression of Ubl3 gene across all tissue types and replicates*

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/Nom1Expression.png)

*Figure 5: Expression of Nom1 gene across all tissue types and replicates*


### Enhancer analysis 

Figure 6 captures some of the highest H3K27Ac peaks in the data set, accompanied by H3K4me1 peaks. Across the 4 genes covered by this stretch that have high expression (Ppm1g, Zfp513, Snx17, Eif2b4), H3K27Ac and H3K4me1 peaks occurred at the beginning of them. Figure 7 shows Cenpa gene region where high expression is accompanied by H3K27Ac and H3K4me1 peaks at the beginning of the gene with H3K4me1 stretching into the introns. In the case of Dpysl5 gene shown in figure 8, mid brain tissue had higher levels of H3K27Ac and H3K4me1 than fore limb and hind limb tissues. Peaks started at the beginning of the gene stretching across intronic regions of the as well. In this case, the mid brain with higher histone modification peaks stretching for further distances than the hind and fore limbs has much greater levels of exon region expression. This highlights the correlation that regions of higher H3K27Ac and H3K4me1 abundance result in greater expression of genes immediately downstream from them, classifying H3K27Ac and H3K4me1 as enhancers. Peaks in the PhyloP track correspond to regions that are highly conserved across different species. For the most part, exonic regions have the highest peaks. In the case of the Shh gene, some intronic regions have peaks as high as the exonic ones. This may correspond to regulatory regions of other genes that are also highly conserved. 

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/fig6VariousGenes.png)

*Figure 6: IGV window showing Ppm1g, Zfp513, Snx17, and Eif2b4 regions with expression and histone modification data in different tissue types*

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/Cenpa.png)

*Figure 7: IGV window showing Cenpa gene region expression and histone modification data in different tissue types*

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/Dpysl5.png)

*Figure 8: IGV window showing Dpysl5 gene region expression and histone modification data in different tissue types*

In IGV chr5:29,314,718-29,315,770 location (seen in figure 9) corresponds to Zone of polarizing activity regulatory sequence (ZRS) also known as MFCS1, a region known to regulate the Shh gene. This enhancer region had H3K4me1 peaks stretching across its area in all tissues. However, in mid brain the peak was very small. In terms of H3K27Ac, there were little to no peaks in mid brain tissue while hind limb and forelimb had slight peaks. Peaks remained at a relatively constant height across the region with no dramatic spikes. There was little histone appearance in mid brain but high Shh expression so these histones and ZRS region might have a smaller regulatory effect on mid brain cells. The moderate peaks in hind limb and fore limb were consistent with the moderate Shh expression in these tissue types. 

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/ZRSregion.png)
*Figure 9: IGV window showing ZRS regulatory region expression and histone modification data in different tissue types*

The alignment of ZRS enhancer regions from human, mouse, cow, dolphin, python, rattlesnake, cobra, and boa were viewed using mafft. The most interesting region of alignment came in the 241 to 400th bases as shown in figure 10. Human sequence `ttttaggtaacttcc` was well conserved in organisms with limbs but absent in snakes. Human sequence ` tgctggtgcttggaatg` while less neatly conserved in organisms with limbs was also missing in snakes ZRS enhancer region. 

![alt text](https://github.com/cse185-sp18/cse185-week4-notmaurox/blob/master/labreport/mafft.png)
*Figure 10: region 241-400 of multiple sequence alignment of ZRS enhancer region from different species viewed in mafft*

### Extra Credit: Gene ontology analysis
Table 3 contains information on the biological processes top differentially expressed genes are involved in. This data was gathered using the PANTHER classification system [9]. Function reported corresponds to listed PANTHER GO-slim molecular function for each gene. 

|  transcript ID        | Gene Name | Molecular Function |
| --------------------- | --------- | ---------- |
| ENSMUST00000075453.8  | Rpl21     | structural constituent of ribosome |
| ENSMUST00000031131.10 | Uchl1 | cysteine-type peptidase activity |
| ENSMUST00000031249.7  | Sparcl1 | calcium ion binding, protein binding |
| ENSMUST00000056355.8  | Nat8l |	acetyltransferase activity |
| ENSMUST00000079324.13 | Ubl3 | ubiquitin-protein ligase activity |
| ENSMUST00000102553.10 | Hmgn2 | chromatin binding, nucleic acid binding |

*Table 3: Top differentially expressed genes identified by sleuth and process they are involved in identified using PANTHER* 

## Discussion 
Based off the results above, H3K27Ac and H3K4me1 histone modifications correspond to increased gene expression. This is true for many genes across different tissue types in the developing mouse, the primary exception being with the Shh gene. ZRS is understood to be a regulatory region for the Shh gene. However, in mid brain tissues, the sample with the highest Shh gene expression, H3K27Ac and H3K4me1 were seen at very low levels in the ZRS region. This points to the conclusion that in the mid brain, a different regulatory region is determining Shh expression. In the fore limb and hind limb tissues, moderate H3K27Ac and H3K4me1 abundance in ZRS region was correlated with moderate Shh expression, confirming that in limb tissues, the ZRS enhances Shh expression. When ZRS region was compared across species, snakes with no limbs had the least preserved sequence with major chunks missing. This potentially means that they have a loss of function mutation that prevents Shh expression to the level it appears in species with limbs.

Upon review of IGV, the region surrounding En2 and Cnpy1 genes has greater H3K4me1 abundance in mid brain tissue and corresponds to greater expression of these genes in those cells. Slightly higher H3K27Ac levels also accompany this region. The inverse is true for the neighboring Rbm33 region where more histone abundance correlates to lower gene expression. In the case of the Rbm33 gene, PhyloP track shows it as more conserved than En2 and Cnpy1. To identify other regions of limb specific enhancers with high conservation across species, one could filter out regions of low conservation (below a threshold set by average height) in PyloP and in the remaining regions do local searches for limb specific histone modifications in that area. A filtered search process would cut down on computation time and cost.  

In limb tissue in the mouse sample, Shh seemed to play an important role. Furthermore, limbless snakes had the least preserved Shh gene, pointing to the idea that the snake copy has mutated and lost some function. To further confirm that Shh mutations in snakes lead to a loss of legs, one could genetically engineer mice with the snake shh gene and study the resulting phenotype. Alternatively, a mouse Shh gene could be engineered into a snake to see if the resulting offspring have legs. 


## Citations
1. Andrews S. [Internet]. FastQC A Quality Control tool for High Throughput Sequence Data. Babraham Bioinformatics; [cited 2018Apr10]. Available from: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ EcoliWiki [Internet].
2. Harold J. Pimentel, Nicolas Bray, Suzette Puente, Páll Melsted and Lior Pachter, Differential analysis of RNA-Seq incorporating quantification uncertainty (sleuth), Nature Methods (2017), advanced access http://dx.doi.org/10.1038/nmeth.4324.
3. James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011)
4. Yamada, Tomii, Katoh 2016 (Bioinformatics 32:3246-3251) additional information Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. (explains some options for aligning a large number of short sequences)
5. N P Brown, C Leroy, C Sander; MView: a web-compatible database search or multiple alignment viewer., Bioinformatics, Volume 14, Issue 4, 1 January 1998, Pages 380–381, https://doi.org/10.1093/bioinformatics/14.4.380
6. N. L., Bary, et al. “Near-Optimal Probabilistic RNA-Seq Quantification.” Nature Biotechnology, vol. 34, 2016, pp. 525–527., doi:doi:10.1038/nbt.3519.
7. Wei, Iris H., et al. “RNA-Seq Accurately Identifies Cancer Biomarker Signatures to Distinguish Tissue of Origin.” Neoplasia, vol. 16, no. 11, 2014, pp. 918–927., doi:10.1016/j.neo.2014.09.007.
8. Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras; STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, 1 January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635
9. PANTHER version 11: expanded annotation data from Gene Ontology and Reactome pathways, and data analysis tool enhancements.
Huaiyu Mi, Xiaosong Huang, Anushya Muruganujan, Haiming Tang, Caitlin Mills, Diane Kang, and Paul D. Thomas 
Nucl. Acids Res. (2016) doi: 10.1093/nar/gkw1138



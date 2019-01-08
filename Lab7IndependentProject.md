# Kakapo Mitogenome Analysis For SNP Identification: Investigating Potential Markers of Diversity Within Modern and Historical Populations 

#### Mauro Chavez
#### A12150388
#### 6/10/2018

## Abstract

The goal of this report was to identify SNP markers within the mitogenome of the endangered Kakapo parrot to propose a new measurement for cataloging diversity. This involved aligning reads from historical and modern populations to a reference genome, identifying SNPs, and performing filtering to identify significant ones. This is important because in small populations suffering from dramatic population loss, it is hard to catalog genetic diversity, a crucial statistic for increasing the fitness of a recovering population. In the end, 64 SNPs were identified as potential indicators separating the modern, less diverse, population from the historic, more diverse, population. However, poor performance of analysis tools and questionable statistical significance prevents these results from being conclusive. As a result, other methods of analysis and other sources of SNP identification are proposed. 

## Intro

The Kakapo, a flightless parrot, was once one of the most abundant bird species in New Zealand as recently as European arrival to the island [7]. However, today the bird species numbers to 150 members [7]. While there are major research and conservation efforts taking place in New Zealand focused on preserving and rescuing the kakapo species, Ian Jamieson points out some troubling philosophies surrounding this work. He reports, “in New Zealand… it has been claimed that native birds are less affected by inbreeding depression… [due to] a long history of small population size.. [that] has “purged” deleterious alleles” [8]. Furthermore, he claims “this has lead to a general neglect of inbreeding as a factor in recovery programs” [8]. Additionally, researchers in this field have acknowledged the recent dramatic population decline and what that means for genetic diversity [6]. Due to this concern, the goal of this analysis was to determine another metric for tracking diversity in the remaining Kakapo population. Two major studies cataloging genetic diversity and piecing together Kakapo ancestral history have been completed by Dussex et al. and Bergner et al. Both these studies utilize historical Kakapo skin samples and the modern population to generate mitochondrial sequencing and microsatellite data. The consistent methodology across these two studies being that this data can help create an understanding of the “timing and magnitude of population bottlenecks”[6] to gain a better sense of the specific selective pressures put on a species.  

Recent studies have shown that mitochondrial DNA can catalogue the accumulating of deleterious mutations in species experiencing inbreeding [7]. It is with this understanding that the goal of Wright et al’s. study on Tasmanian devils [9] is applied to the mitochondrial sequencing data set provided by Dussex et al. The Tasmanian devil is another endangered species where cataloging diversity is extremely important. While the population decline of the Kakapo has been attributed to European arrival to New Zealand [6], Tasmanian devils have suffered a sharp population decline due to a contagious tumor disease [9]. Wright et al pointed out that when genetic diversity is low in a captive population, “a more sensitive genotyping assay [than microsatellite markers] is required" [9]. Their solution was to generate a profile of unique SNPs within the population to assess genetic diversity. Due to the success of their approach, in this analysis I hoped to accomplish something similar, identify SNPs within the mitochondrial genomes of Kakapo’s that serve as markers of biodiversity. The analytical approach, similar to that of a GWAS study, was focused on creating a profile of mitochondrial SNPs that differentiate the more diverse historical samples from the less diverse modern samples.

## Methods
### Data Set

The data set used in this analysis was collected by Dussex et al for their paper “Full Mitogenomes in the Critically Endangered Kākāpō Reveal Major Post-Glacial and Anthropogenic Effects on Neutral Genetic Diversity” [7]. They used 39 historical Kakapo skin samples from various museum collections to establish a past kakapo population.  79 blood samples from the current protected Kakapo population were used to establish the modern group [7]. The European Nucleotide Archive’s FTP downloader [10] was used download data from 122 samples studied by Dussex et al. 40 historical sample’s mitochondrial DNA was provided as single read libraries in fastq file format.  82 modern sample’s mitochondrial DNA was provided as both single read and paired end read libraries. Dussex et al. used SeqPrep 1.1 to trim adapters and merge paired-end reads to generate single read fastq libraries for the modern samples. For simplicity, only single read libraries were used for this analysis. While Dussex et al. reported a total of 118 samples used for analysis, 122 samples were uploaded to The European Nucleotide archive. This discrepancy seems to be a result of some modern samples being sequenced twice. For example, Kakapo Sarah appears as samples `Sarah` and `Sarah_new`. For the sake of a more diverse data set, all samples uploaded to The European Nucleotide archive were used in analysis. 

### Analysis Procedure 

The majority of the analysis pipeline was automated by **analysisScript.sh.** In order for this script to work, a csv containing sample IDs and population labels was created to distinguish historical (H) from modern (M) samples. The first step of this script was to build a list of all fastq files and concatenate them into one file. This file was used by **FastQC** v0.11.07 [5] in order to gain statistics on read qualities across all 122 samples.  Once a FastQC report was generated, the script proceeded to align all fastq files to reference genome NC_005931.1. **BWA** version 0.7.12-r1039 [2] was used to index the reference genome prior to use of analysisScript.sh. For each sample’s single read file, bwa mem was used to align reads to reference and create .sam files. From there, **SAMtools** v 1.5 [3] was used for compression into .bam files and then for sorting and indexing resulting bam files. Above steps were completed with default parameters. Once .bam files were created for all samples, they were used by SAMtools to create an mpileup file. Additionally, mpileup files for historic and modern Kakapo populations were created as well. `analysisScript.sh` covered automation of all the above steps. 

**VarScan** version 2.3 mpileup2snp [4] was used outside analysisScript.sh to calculate the number of rare variants ( `--min-var-freq 0.001` ) and common variants (  `--min-var-freq 0.95` ) in the entire data set but also the historic and modern subsamples. Mpileup2snp was run with ` --variants ` and ` --output-vcf 1 ` flags. 

Once .vcf files were created, **PLINK** version 1.90b5.4 was used to process the SNPs identified as common variants in the entire data set. The motivation being to filter out SNPs that were common across the historic and modern samples that could be used to separate the two groups. ` plink --vcf ` was used to process .vcf files into .bed, .bim, and .fam files with `--allow-extra-chr --out ` flags.`--allow-extra-chr` flag allowed for chromosome number to be represented by reference genome ID. The same tool was used with `--maf 0.05 –allow-extra-chr --recode` flags to make .ped and .map files. A .phen file was created manually by altering the data present in the .csv file used in analysisScript.sh. Due to the relatively low number of SNPs that came out of VarScan, principal component analysis was skipped in order to preserve more SNPs that could potentially serve as measurements of diversity. Plink logistic regression was completed with flags `--allow-no-sex`, `--allow-extra-chr`, and `--pheno` to specify the phenotype file. This resulted in an `assoc.logistic` file that contained the final list of potentially significant SNPs in tracking diversity in historical vs modern kakapo populations. Additionally, a linear regression model was tested as well with `--all-pheno` and `--linear` flags. 

### analysisScript.sh code
```
#!/bin/bash
#From list of all samples, concatinate reads into single fastq file for fastqc
INPUT=KakapoPopList.csv
OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
totalSampleList=""
while read sample type
do
 totalSampleList+="KakapoExpanded/$sample.fastq "
done < $INPUT
echo $totalSampleList
IFS=$OLDIFS
cat $totalSampleList > allSamples.fastq
fastqc -o . allSamples.fastq
# Perform alignment and generate sorted bam files for all samples. 
END=2519160
outPutString = ""
for((x = 2519050;x<=END;x++));
do
 echo "starting  ERR$x"
 bwa mem NC_005931.1.fasta KakapoExpanded/ERR$x.fastq > ERR$x.sam
 samtools view -S -b ERR$x.sam > ERR$x.bam
 samtools sort ERR$x.bam > sERR$x.bam
 samtools index sERR$x.bam
 outPutString+="sERR$x.bam "
 echo "finished ERR$x"
done
for y in 2515104 2515103 2510938 2510937 2510936 2510935 2510934 2506007 2506006 2505730 2505729
do
 echo "starting  ERR$y"
 bwa mem NC_005931.1.fasta KakapoExpanded/ERR$y.fastq > ERR$y.sam
 samtools view -S -b ERR$y.sam > ERR$y.bam
 samtools sort ERR$y.bam > sERR$y.bam
 samtools index sERR$y.bam
 outPutString+="sERR$y.bam "
 echo "finished ERR$y"
done
#Make an mpile up for full data set
samtools mpileup -f NC_005931.1.fasta $outPutString > combKakapo.mpileup
#Create mpileup for historic and modern samples
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
historicList=""
modernList=""
while read sample type
do
 if [[ $type == *"M"* ]]; then
  modernList+="s$sample.bam "
 fi
 if [[ $type == *"H"* ]]; then
  historicList+="s$sample.bam "
 fi
done < $INPUT
IFS=$OLDIFS
samtools mpileup -f NC_005931.1.fasta $historicList > historicSample.mpileup
samtools mpileup -f NC_005931.1.fasta $modernList > modernSample.mpileup
```

## Results

### FastQC Results

FastQC was used to compute quality control metrics on a file containing all single sequencing reads across all samples, some of this data is reported in Figure1. Single reads from modern samples represent paired end reads that had adaptors trimmed and reads merged by Dussex et al with SeqPrep 1.1 [7]. Per tile sequence quality module failed where red tiles only appeared at the very end of reads in the 235-237 bp range.  Apart from that, the graph was overwhelmingly blue. The sequence length distribution module threw a warning. The graph appears normal with a peak at the 50 bp sequence length followed by a steady decline, but at the ~130 bp marker there is another sharp dramatic spike that appears to go higher than the first peak.  The sequence duplication level module threw a warning because non-unique sequences make up more than 20% of total reads. Other than that all modules passed FastQC analysis. 

![alt text](https://github.com/cse185-sp18/cse185-final-project-notmaurox/blob/master/Screen%20Shot%202018-06-10%20at%205.05.37%20PM.png)

*Figure1: FastQC analysis of all reads across all samples* 

### VarScan Variant Identification

Analysis script `analysisScript.sh` created mpileup files for 3 groups, the complete data set, the historical kakapo samples subset, and modern kakapo samples subset. VarScan was used to identify common variants, with minimum variant frequency 0.95, and rare variants, with minimum variant frequency 0.001, in all three groups. Data for common variant identification is reported in Table1 and data for rare variant identification is reported in Table2. The Complete data set and historic samples seemed to have coverage across the entire reference genome while modern samples were missing ~30 bases. In the historic samples, the strand filter did not remove any variants. As to be expected, when searching for rare variants, there were more positions reported than if searching for common variants. Interestingly, relative to the size of the reference genome, SNPs made up slightly less than 1% of the genome. 

|   |  Number of  bases in pileup file | Common Variant Positions | Common Variant Positions Failed by strand-filter | Reported Common Variants|
|---|---|---|---|---|
| Complete Data Set | 15599 | 220<br>(218 SNP, 2 indel) | 39 | 179<br>(179 SNP) |
| Historic Samples | 15599 | 100<br>(99 SNP, 1 indel) | 0 | 99<br>(99 SNP) |
| Modern Samples | 15562 | 126<br>(125 SNP, 1 indel) | 41 | 84<br>(84 SNP) | 

*Table1: VarScan statistics for common variant identification `--min-var-freq = 0.95`*

|   |  Number of  bases in pileup file | Rare Variant Positions | Rare Variant Positions Failed by strand-filter | Reported Rare Variants|
|---|---|---|---|---|
| Complete Data Set | 15599 | 286<br>(280 SNP, 7 indel) | 41 | 239<br>(239 SNP, 1 indel) |
| Historic Samples | 15599 | 127<br>(123 SNP, 5 indel) | 0 | 123<br>(123 SNP, 1 indel) |
| Modern Samples | 15562 | 167<br>(165 SNP, 2 indel) | 43 | 122<br>(122 SNP) | 


*Table2: VarScan statistics for rare variant identification `--min-var-freq = 0.001`*

### PLINK Results

The goal of PLINK analysis was to identify the significant SNPs that differentiate the historical group from the modern one. The first step was to convert the .vcf file into the files required for PLINK. The statistics reported in Table3 capture the modifications made to the common variant data across all samples in this process. Minor allele frequency cut off was set to 0.05 (5%) to give greater confidence that these variants were prevalent in the Kakapo communities. This cut off resulted in 94 variants being removed from the list of SNPs. The genotyping rate of 0.410202 reveals that only 41% of SNPs were useful in determining genotypes of samples.  In the case of most SNPs, across other samples, they were missing. In the end, 85 variants remained. 

The resulting files were used for both logistic and linear regression. Quantitative traits are recommended to be studied using linear regression while disease traits are recommended to be studied using logistic regression. Both models were tested to see how the results differed. `SupplementaryData.md` contains supplementary tables 1 and 2 that report statistics on the SNPs identified by both models. In both cases, the full set of 85 SNPs were reported in the results. However, in the case of the logistic model (supplementary table 1), 21 SNPs had no odds ratio, coefficient t-statistic, or asymptotic p-value for t-statistic reported. The same is true for results from the linear model (supplementary table 2). In both models, the historic samples served as the controls and the modern samples served as the cases, and in both models the same results were reported. 

| | Variants loaded from .bim file | Num samples | Genotyping rate | Variants removed due to minor allele thresholds | Reported Variants |
|-|--------------------------------|-------------|-----------------|-------------------------------------------------|-------------------|
| Complete Data Set Common Variants | 179 | 122 | 0.410202 | 94 | 85 |

*Table3: Results from processing common variant .vcf file into files required for PLINK with --maf 0.05 parameter*

## Discussion

The results of this report leave mitochondrial SNPs in a questionable state when it comes to their usefulness as markers of genetic diversity identified using the approach followed in this report. FastQC data was consistent with knowledge of how the single end reads were generated. The majority of them were created by merging paired end reads, which explains why quality scores remained so high near the end of the reads. 40 historical samples were sequenced as single end reads which explains why the read length curve appears the way it does. The second peak must represent the artificial single reads that were the result of merging paired end reads. 

VarScan data reported more variants when scanning for rare ones as opposed to common ones. However, relative to the amount of common variants identified, not many more appeared once the minimum variant frequency was dramatically lowered. This might represent a biological constraint of mitochondrial genomes where they are unable to acquire many SNPs. Data on variants identified is also influenced by the way in which the reference genome was assembled. In the case of the endangered Kakapo, there arent many samples available for the creation of a reference genome that better represents the species as a whole. On average, modern samples had more SNPs than historic samples and this bias might be the result of the data that went into creating the Kakapo mitochondrial reference genome. 

Further statistical analysis of results was skipped due to unfamiliarity with specifics of mitochondrial DNA behavior in addition to small set of SNPs provided as results. The goal was to complete a preliminary investigation of potential SNP markers that could be used in diversity assays for Kakapos, so as many as possible were kept in the running at all steps. The statistics reported during .vcf conversion into files necessary for PLINK computations, point to potential issues with this data set. The low genotyping rate and dramatic reduction of SNPs that remained after conversion reveal the inability of these SNPs to present a holistic picture of the diversity within the data set. When it comes to genotyping by SNPs, one hopes for many more to be identified with a higher genotyping rate to increase the chance SNPs are found that effectively distinguish between groups of study. Additionally, in the final logistic and linear regression files, data was missing from OR/STAT/P fields that would be required for statistical analysis. The fact that the results from both the linear and logistic models were identical raises suspicions that this pipeline failed to capture meaningful data. Principal component analysis was attempted multiple times but in each case, all zero eigen vectors were returned. This is representative of problematic aspects of the pipeline used in this report. For further investigation of this issue, an alternate analysis pipeline could be designed using the BAM files of complete mitogenomes provided by Dussex et al [7]. However, at the time of completing this analysis, downloading these files from The European Nucleotide Archive resulted in an error.

In the case of Wright et al[9], SNPs were searched for in data resulting from whole genome sequencing. This increase in the amount of data available for analysis greatly contributed to their success in identifying SNPs useful as genetic markers of diversity. However, mitochondrial genomes are much smaller and result in less data to work with. While mitochondrial genomes would be nice to use for analysis due to their already inherit utility in other studies such as Dussex et al[7] and Bergner et al[6], the failure of PLINK to provide meaningful results is representative of the data’s limitations.  To verify these limitations, I propose using the complete mito-genomes provided by Dussex et al and a different analysis pipeline to generate the files required for PLINK analysis. If similar problematic results are obtained, mitochondrial DNA may be classified as not useful for generating SNP profiles for measuring genetic diversity. If the same set of 64 SNPs is returned, it may be worthwhile to reopen the case for mitogenome SNPs as a catalogue of genetic diversity. These results would need to be backed up by building a classifier and testing it against other data sets (hopefully future populations of Kakapo parrots). Ultimately, this procedure followed a GWAS approach to identifying genetic markers distinguishing historical and modern Kakapo populations. This approach may not be appropriate for data sets like the one analyzed in this report. Alternatively, the tools used might signify the need for a more streamlined pipeline. 

## Citations
[1] Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81

[2] Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform (BWA-MEM). Bioinformatics, 25:1754-60. [PMID: 19451168]

[3] Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

[4] Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111

[5] Andrews S. [Internet]. FastQC A Quality Control tool for High Throughput Sequence Data. Babraham Bioinformatics; [cited 2018Apr10]. Available from: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ EcoliWiki [Internet].

[6] Bergner LM, Dussex N, Jamieson IG, Robertson BC. European Colonization, Not Polynesian Arrival, Impacted Population Size and Genetic Diversity in the Critically Endangered New Zealand Kākāpō. Journal of Heredity. 2016Feb;107(7):593–602.

[7] Dussex N, Seth JV, Robertson B, Dalén L. Full Mitogenomes in the Critically Endangered Kākāpō Reveal Major Post-Glacial and Anthropogenic Effects on Neutral Genetic Diversity. Genes. 2018;9(4):220.

[8] Jamieson IG, Wallis GP, Briskie JV. Inbreeding and Endangered Species Management: Is New Zealand Out of Step with the Rest of the World? Conservation Biology. 2006Sep;20(1):38–47.

[9] Wright B, Morris K, Grueber CE, Willet CE, Gooley R, Hogg CJ, et al. Development of a SNP-based assay for measuring genetic diversity in the Tasmanian devil insurance population. BMC Genomics. 2015;16(1).

[10] Downloading read and analysis data [Internet]. ENA FTP Downloader. European Nucleotide Archive; Available from: https://www.ebi.ac.uk/ena/browse/read-download#downloading_files_ftp



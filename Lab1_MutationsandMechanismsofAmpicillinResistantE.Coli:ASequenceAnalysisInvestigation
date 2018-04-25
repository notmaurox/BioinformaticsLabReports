# Mutations and Mechanisms of Ampicillin Resistant E. Coli: A Sequence Analysis Investigation  
##### Mauro Chavez

## Abstract
Analysis of next generation sequencing data of ampicillin resistant E. coli strain provided a detailed report of mutated genes potentially causing this behavior. A fastqc, sickle, BWA-MEM, samtools, VarScan, VEP workflow was used to enhance the quality of raw paired end reads before filtering so that alignment to a reference genome and analysis of mismatched bases revealed information on bacterial cell mutations. Mutations were searched against databases to discover the altered mechanisms of the resistant E. coli strain, a process critical to a world where the gaining of antibiotic resistance threatens major advances in medicine extending human life expectancy[8]. 

## Introduction
This project involved the analysis of an ampicillin resistant MG1655 E.coli strain through studying sequencing data, a workflow that has great potential benefits in clinical settings where doctors are faced with patients suffering from bacterial infection. In her article “Health care: Bring microbial sequencing to hospitals”, Sharon Peacock touches on the current state of personalized drug combinations to combat drug resistant bacteria.  It would take up to 8 weeks of lab work to determine a plan of attack while a doctor makes an educated guess as to what might be helpful until lab results are returned [8]. Alternatively, it could be less than a week before a doctor is provided sequence data for a culture of bacteria from their patient [8].  With this data, Peacock imagines a doctor could use specially designed automated software to identify mutations found in the sample from a global microbial genome database and determine the optimal treatment plan. 

This is not only an important tool for combating for current drug resistant bacteria, but also for fighting future drug resistant bacteria. Peacock reminds us that in the 1930’s a baby had a life expectancy of 60 years, a number that has grown to 80 years today. However, she warns that antimicrobial resistance threatens the progresses made in human longevity aided by antibiotics [8]. For this reason, it is critical that the analysis done in this report becomes more streamlined and common place in the medical field.
 
## Methods
Analysis began by collecting three files, amp_res_1.fastq and amp_res_2.fastq, and a reference MG1655 E. Coli genome NC_000913.3.fasta. These fastq files contained foreword and reverse paired end read sequencing data from an ampicillin resistant E. coli strain. Fastqc[6] was run on these two files with default parameters. From the output html files provided by fastqc analysis, the “Per base sequence quality” and “Per tile sequence quality modules” failed. This meant that the bases with lowest quality scores exceeded a threshold defined by the quality of all bases in the same data set. 

To filter the reads, Sickle[1] was used to trim them within a certain quality and length threshold. Sickle was run with default parameters and its output trimpair1.fastq and trimpair2.fastq were reassessed by fastqc. The trimmed reads no longer failed the quality score modules within the fastqc html files. 

The trimmed reads (foreword and reverse) were aligned to reference genome NC_000913.3.fasta using BWA-MEM[2] and the output was redirected to a .sam file. This sam file was decoded by flagstat, a samtool utility [3]. Samtools was then utilized to compress and sort the sam file, creating an indexed bam file. This bam file was viewed using tview and mapped reads were compressed using mpileup, both utilities of samtools. The compressed file provided by mpileup was then used by VarScan[4] with a min var freq of 0.50 to create a .vcf file. An awk script as defined in the lab report and notebook helped trim unnecessary data and format columns of this .vcf file so that it could be uploaded to the online Variant Effect Predictor [5] (E. coli MG1655) to locate mutations in the genome and discover which genes they fell in. The awk script removed the first 24 lines and replaced the first cell of each row with keyword “Chromosome” and outputted the edited file to  a new one. VEP data was viewed online in the web browser used for submitting the processed .vcf file. 

## Results

The fastqc per base sequence quality analysis of amp_res_1.fastq is presented in Fig 1. The per base sequence quality for amp_res_2.fastq very closely resembled that of amp_res_1.fastq except on average it had wider error bars. 

![alt text](https://github.com/cse185-sp18/cse185-week1-notmaurox/blob/master/labreport/CSE185Lab1PreTrim.jpg)

*Fig 1: Per-position read quality before trimming by sickle of amp_res_1.fastq*

After trimming using Sickle, quality scores increased dramatically. Fig 2 contains information on how the reads were trimmed. 

```
PE forward file: /home/linux/ieng6/cs185s/public/week1/amp_res_1.fastq
PE reverse file: /home/linux/ieng6/cs185s/public/week1/amp_res_2.fastq

Total input FastQ records: 7107040 (3553520 pairs)

FastQ paired records kept: 6904494 (3452247 pairs)
FastQ single records kept: 98694 (from PE1: 94870, from PE2: 3824)
FastQ paired records discarded: 5158 (2579 pairs)
FastQ single records discarded: 98694 (from PE1: 3824, from PE2: 94870)
```

*Fig 2: Sickle output describing how foreward and reverse reads were trimmed*

Fig 3 contains the resulting per base sequence quality scores after trimming of amp_res_1.fastq. All the bases now lie within the green range with the exception of the last nucleotide. Similar to before trimming, amp_res_2.fastq after trimming by sickle looked very similar to the quality scores of amp_res_1.fastq except with wider error bars. In both cases, the error bars for the last nucleotide provide the potential that its score falls into the red zone. However, given that these are paired end reads, the last nucleotide of the foreword read might correspond to the first nucleotide of the reverse read. 

![alt text](https://github.com/cse185-sp18/cse185-week1-notmaurox/blob/master/labreport/CSE185Lab1PostTrim.jpg)

*Fig 3: Per-position read quality after trimming by sickle of amp_res_1.fastq*

Table 1 corresponds to the number of reads before and after trimming. Additionally, it reports how many reads were aligned to the reference genome, data that was provided by using flagstat utility of samtools on the .sam file that reported the alignments of the trimmed reads to the reference NC_000913.3.fasta genome. 

| Initial num. reads | Num reads after trimming | Num. reads aligned |
| -------------------| ------------------------ | ------------------ |
| 7107040            | 6904494                  | 6899109            |
| Initial num. paired end reads   |  Num. paired end reads after trimming           |                    |
| 3553520 | 3452247 |                    |

*Table 1: Number of reads and paired reads before and after trimming, including final number of reads aligned*

Once the .sam file was compressed and sorted into a bam file, it was processed into an mpileup file. VarScan was used to single out variants in the piled reads and outputted into a .vcf file. After some processing, this file was given to the Ensemble Variant Effect Predictor. The data from VEP processing is reproduced in Table 2 followed by descriptions of the genes effected by missense mutations. 


| Location         | Consequence   | In Gene?| Gene Name | AA subsitiution | Codeon Change   |
| ---------------- |:-------------:| -------:| ---------:| ---------------:| ---------------:|
| 92439-92439      | missense      | Yes     |    b0084  |    D/N          | **G**AC/**A**AC |
| 803662-803662    | missense      | Yes     |    b0771  |    L/I          | **C**TT/**A**TT |
| 852762-852762    | -             | -       |    -      |    -            | -               |
| 1905761-1905761  | missense      | Yes     |    b1821  |    G/D          | G**G**T/G**A**T |
| 3535147-3535147  | missense      | Yes     |    b3404  |    V/G          | G**T**A/G**G**A |
| 4390754-4390754  | synonymous    | Yes     |    b4161  |    A/A          | GC**C**/GC**A** |

*Table 2: http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Tools/VEP [5] identified mutations*

**Gene:**  b0084 | ftsl

Function: Penicillin-binding protein 3, also functions as transpeptidase aiding in septal peptidoglycan synthesis [5]. Essential for cell division described by ecoliwiki.net as “transpeptidase involved in septal peptidoglycan synthesis (penicillin-binding protein 3)” [7].

**Gene:** b0071 | ybhJ  

Function: Predicted hyrdatase [7] of aconitase family of proteins [5].

**Gene:** b1821 | mntP

Function: Described by EnsemblBacteria as “putative Mn(2+) efflux pump, mntR-regulated” [5] and does not appear in ecoliwiki.net database. 

**Gene:** b3404| envZ

Function: Both EnsemblBacteria and ecoliwiki.net agree that this is a “sensory histidine kinase in two-component regulatory system with OmpR” [5][7].

## Discussion
It is expected that as read length progressed that quality scores decreased due to technical errors of sequencing procedures. Luckily having paired end reads allowed for the making up of this decrease in quality in the last nucleotide of each read. The reads themselves were greatly improved by sickle.

Over 99% of the trimmed reads were aligned showing the effectiveness of sickle to create useful map-able reads. However, in the data provided by VEP, there were some unusual results. For example, the third nucleotide variant provided no data on the mutation type. This is potentially a result of the very loose N parameter used in the VarScan step with the mpileup file in creating the .vcf file. N = 0.5 meant that at least 50% of non reference bases at a position did not have to match the reference for classification as a mutation. Given a larger N, this third reported variant would most likely disappear. A similar story might be true or the synonymous mutation identified as the last mutation. 

The first reported missense mutation in b0084 was the most interesting. As a gene involved in the production of a penicillin binding protein, it helps define the bacteria’s interaction with the drug.  The mutation in this gene could potentially sever communication with penicillin allowing the bacteria to ignore its presence. The third missense mutation occurred in gene b1821, transcribing a manganese pump. This mutation could potentially lead to a gain of function mutation allowing the bacteria to better defend against drugs that block secretion of manganese ions. The forth missense mutation occurred in b3404, a gene that codes for a kinase that is apart of a regulatory system. UniProt describes this gene as coding for a protein that regulates osmoregulation and phosphorylates in response to environmental signals[9]. This may contribute to a bacteria’s ability to better respond to threats posed by the presence of specific drugs or ignore fake signals given off by drugs tricking the bacteria into shutting down. 

I would recommend a treatment for this E. coli that avoids penicillin related attack mechanisms. Additionally, I would recommend a drug that’s target is not a manganese related pathway because the one functioning in this E. coli may differ dramatically from that of other types of E.coli. Ultimately, sequencing provided a glimpse into the specific cellular mechanisms that are being disrupted in the ampicillin resistant bacteria, information that would otherwise be acquired from much more time intensive and rigorous testing of the mutated bacteria. 

## Citations
1. Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files (Version 1.33) [Software].  Available at https://github.com/najoshi/sickle.
2. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
3. Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
4. Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111 
5. McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4 Available at http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Tools/VEP
6. Andrews S. [Internet]. FastQC A Quality Control tool for High Throughput Sequence Data. Babraham Bioinformatics; [cited 2018Apr10]. Available from: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ EcoliWiki [Internet]. 
7. EcoliWiki:About - EcoliWiki. PortEco; [cited 2018Apr10]. Available from: http://ecoliwiki.net/colipedia/index.php/EcoliWiki:About 
8. Peacock S. Health care: Bring microbial sequencing to hospitals. Nature [Internet]. 2014;509(7502):557–9. Available from: https://www.nature.com/news/health-care-bring-microbial-sequencing-to-hospitals-1.15282 
9. European Bioinformatics InstituteProtein Information ResourceSIB Swiss Institute of Bioinformatics. "Osmolarity Sensor Protein EnvZ." EnvZ - Osmolarity Sensor Protein EnvZ - Escherichia Coli (strain K12) - EnvZ Gene & Protein. March 28, 2018. Accessed April 10, 2018. http://www.uniprot.org/uniprot/P0AEJ4.


# De Novo Assembly Pipeline: Building Staphylococcus Genome from Short Reads
##### Mauro Chavez

## Abstract

This week in lab we worked through a de novo assembly pipeline to piece together a Staphylococcus aureus genome using two libraries of paired end reads. This species is important not only because it is commonly used to benchmark assembly protocols but also because it is involved in antibiotic resistant bacteria. Reads were assembled into contigs by identifying regions of overlap. Paired end reads with large insert size were then mapped to separate contigs to arrange them into scaffolds. Reads were then remapped to these scaffolds to fill in gaps ultimately resulting in a 97.865% assembled genome.

## Introduction

De novo assembly is the process of using short sequencing reads to build a reference genome. It is needed for studying the genomes for species that do not have previously established working reference genomes. Researchers run into this when they study non-model organisms, survey new species when doing conservation research, or survey divergent species. Being able to assemble a genome is also important when studying genomes that have large structural variations from their published reference. Reads from large insertions of nucleotides not present in a reference would need to be built into the reference using de novo assembly. Furthermore, when there are unculturable species within a host, sequencing a sample, sorting out host DNA, and assembling the remaining reads can be used to uncover the genome of said species. 

Assembly requires four main steps. The first is contig assembly where short reads with regions of overlap are conjoined into longer stretches of genome regions. From there, contigs are combined with paired end reads with large insert sizes between them to gain a sense of how contigs are oriented in the genome. Remapping reads to the regions of the genome assembled thus far filled gaps between scaffolded regions. In the case that a reference genome already exists, the assembly can be mapped to the reference to gain a sense of how good the assembly was. 

The goal of this lab was to assemble a Staphylococcus aureus genome from two libraries of short reads. The first library, denoted with “frag” prefix, contained paired end reads with reads oriented towards each other. These reads were 101 bp long and 180 bp apart. These “frag” reads were used to assemble contigs. The second library, denoted with “short_jump” prefix, contained short 37bp reads oriented away from one another with a gap of 3500bp between them. These were used for scaffolding. 

In their paper “Sequence assembly demystified”, Nagarajan et al describe some of the challenges involved in de novo assembly. They explain “Most mathematical formulations of assembly suggest that finding the optimal assembly could require prohibitive computational resources…optimal solutions may be possible for some assembly tasks, such as scaffolding and finishing; however, substantial work still remains to be done before efficient and optimal ‘black box’ solutions for assembly are available” [10]. Here they touch on a critical component to this lab. Unlike previous labs, this week the class is not out to solve a question or entertain a hypothesis. Instead, we all tested different parameters to see who could create the best assembly.  Due to our limited time and computing resources, we could not stress over an optimal assembly, but instead play around with the tools at each step to determine the best approximation. 


## Methods

Once the four fastq files (frag_1.fastq, frag_2.fastq, short_jump_1.fastq, and short_jump_2.fastq) were collected, total number of reads in each file was gathered by counting lines with unix command wc –l and dividing by 4. From there, each file was analyzed using FASTQC [1] to collect quality score data. Quality scores were very poor across all four files. Low quality scores were addressed in files frag_1.fastq and frag_2.fastq using a k-mer based approach. This called for a scan of k-mer distribution in these files using Jellyfish’s count command [9]. Jellyfish counted frequency of all possibly k-mers of provided input length within the dataset. Jellyfish count module was run with –m flag set to 31, counting occurrences of k-mers with length 31. –c flag was included so that directionality of reads was ignored. –s flag was included with estimate genome size of 10000000. The output of this command was used for Jellyfish’s histo command to build a histogram file of k-mer occurrences.  This histogram file, with the help of IPython packages Matplotlib [2] and Pandas[3], was turned into a k-mer distribution graph and saved as a pdf. This graph was used to identify the initial valley point of k-mer distribution, a number that would serve as an input parameter in later steps.

Low frequency k-mers were regarded as likely errors. To get rid of them, SOAPdenovo2[4] was used to scan for similar k-mers (few nucleotide differences) with higher frequency to be used in their correction.  KmerFreq_HA module was run with maximum k-mer size (27), -L flag set to 101 (maximum read length in data set) and –i flag set to 10000000 for both frag_1.fastq and frag_2.fastq. This provided a hash-table that contained the frequency of every k-mer in the data. Corrector_HA module was used with this hash table file with –k set to 27, -Q set to 33, -o set to 3 and –l set to 5, the previously identified valley point in k=31 histogram. Corrected fastq files were outputted with filetype .fastq.cor.pair_[1/2].fq. Jellyfish and IPython steps were repeated with the corrected frag_1 reads to generate a histogram of corrected k-mers. New valley point appeared at 2 and the first peak occurred at 14. 

awkCommand1.txt was used to calculate the average read length. awkCommand2.txt was used to calculate the total number of bases in all the reads. The genome size was then estimated using the previously identified k-mer peak, k-mer size, average read length, and total bases. 

Assembly into contigs was completed using Minia [8] on the corrected read files with k mer size of 29. The class ran minia with kmer sizes ranging from 27 to 43 and resulting statistics were compared at a later point. Repeated Jellyfish analysis showed that first peak remained at 2 with kmer size set to 29 so –abundance-min flag was set to 3 because k-mers with abundance higher than those in valley point were to be used in assembly.  Output files with suffix .contigs.fa contained assembled contigs. Counting lines in the output file and dividing by 2 generated the number of contigs.  The longest and shortest length contigs were recovered using the two commands listed in awkCommand3.txt. QUAST [5] was used to gain some statistics on the contig assembly, skipping contigs shorter than 100, searching for prokaryotic genes, and leaving genome as unknown.

Contig assembly results from class data using kmer sizes 27 to 43, were used to make two graphs, k-mer size vs maximum contig length and k-mer size vs contig N50. 

Reads from short_jump_1 and short_jump_2 were used to scaffold the contigs produced by Mania. For each pair of reads, if one read mapped to one contig and another read mapped to a separate contig, using information about distance between read pairs, distance between contigs and their orientation could be inferred. SSPACE[6] was used for scaffolding. It’s library file was set up with the two short_jump files, insert size of 3500, insert error of 0.5, and orientation set as RF. These parameters aligned with our data being that reads were 3500 base pairs apart and pairs were oriented in reverse – forward direction. SSPACE_Basic_v2.0.pl was run with minimum size required for contigs set to 100 with –s set to the mania contigs file, –p 1 to signify dot file creation and -v 1 to set verbose output. The resulting summary file was studied and statistics were recorded. 

Short_jump read files were cleaned up using Sickle[7]. Sickle was run with –t parameter set to sanger and -l parameter set to 25 on both short_jump files to create two trimmed read fastq files. SSPACE was run again using trimmed short_jump reads for scaffolding and the summary file was again studied for key statistics. 

From the scaffolds built from untrimmed short_jump reads, SOAPdenovo’s GapCloser[4] was used to align reads to scaffolds in an attempt to close gaps between contigs. GapCloser was run with –l file set to 101 and –b set to soap.config. See soap.config file in lab report for full library set up. Resulting assembly was analyzed using QUAST. Contigs were aligned to Staphlyococcus aureus strain NC_010079 reference genome to evaluate de novo assembly performance. Both .fasta and .gff3 from the reference genome were used in QUAST analysis. QUAST was run with scaffolds, find genes, and prokaryotic options checked. A second analysis was run without scaffolds button enabled so that “Icarus”[11] browser could be used in studying assembly. QUAST output was studied to gauge resulting assembly quality.  

## Results

The results of initial FASTQC analysis of the four fastq files included in our reads libraries is displayed in Figures 1-4. Quality scores for the reads in all cases were very low with large error bars. 

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/frag_1_fastqcPerBaseSeqQuality.png)

*Figure 1: frag_1 per base sequence quality*

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/frag_2_fastqcPerBaseSeqQuality.png)

*Figure 2: frag_2 per base sequence quality*

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/short_jump_1_fastqcPerBaseSeqQuality.png)

*Figure 3: short_jump_1 per base sequence quality*

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/short_jump_2_fastqcPerBaseSeqQuality.png)

*Figure 4: short_jump_2 per base sequence quality*

Jellyfish count and histogram analysis provided a .histo file that was used in IPython with the help of Matplotlib and Pandas to create the histogram shown in Figure 5. This was before any read correction had occurred. 

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/frag_1_31kmerHisto.png)

*Figure 5: 31-mer count analysis histogram of frag_1 data. Num kmers vs count*

After SOAPdenovo2 read correction using KmerFreq_HA and Corrector_HA, the resulting .cor.pair_1.fq fastq file containing corrected reads were put through the same Jellyfish IPython pipeline to create the histogram seen in Figure 6. 

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/frag_1_31kmerHistoCorrected.png)

*Figure 6: 31-mer count analysis histogram of frag_1 data after correction. Num kmers vs count*

Table 1 contains the data of the frag fastq files after read correction had occurred. This data was generated using awkCommand1.txt and awkCommand2.txt as mentioned in the method section. 

| corrected file | average read length | total number of bases |
| - | ------------------------------------- | --------------------------------------- |
| frag_1 | 95.6211 | 100818279 |
| frag_2 | 93.3242 | 98396570 |

*Table 1: Shows average read length and total number of remaining bases after awkCommand1 and awkCommand2 usage on correct read files*

Figures 7 and 8 show scatterplots of class data from minia contig assembly with varying kmer size parameters. My data is marked by a red square while class data appears as blue diamonds. Darker diamonds on the plot signify more data supporting that point. 

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/kmerLengthVsMaxContigLengthGraph.png)

*Figure 7: class data for k-mer length vs max contig length, personal data symbolized by square pt*

![alt text](https://github.com/cse185-sp18/cse185-week3-notmaurox/blob/master/labreport/kmerSizeVsContigN50Graph.png)

*Figure 8: class data for k-mer length vs contig N50 score, personal data symbolized by square pt*

The first row of table 2 contains statistics on the assembly of contigs using the corrected reads by minia. The second row contains statistics on the assembly of scaffolds using SSPACE to add short_jump read data to the contigs generated by minia. The third row shows the statistics of scaffolding that used sickle trimmed short_jump reads instead of the raw short_jump read data. The fourth row contains statistics of the “final” alignment that was completed when gaps between scaffolds were filled by alinging of read data to constructed genome using GapCloser. 

| analysis step | number of contigs | max contig length | N50 |
| ------------- | ----------------- | ----------------- | --- |
| Corrected kmer list minia assembly | 740 | 80332 | 25117 |
| After scaffolding | 307 | 232075 | 110285 |
| After scaffolding with sickle output | 516 | 80332 | 25117 |
| After gap filling with untrimmed reads | 307 | 232485 | 110286	|

*Table 2: Contains data on the number of contigs, max contig length, and N50 score at various steps in the assembly proceedure*

Final QUAST analysis provided two outputs, GapCloseroutput and GapCloseroutput_broken. The broken assembly is created by QUAST as a sort of alternate assembly. It uses continuous regions of N’s with length greater than or equal to 10 to split the input fasta into pieces. It then tries to reconstruct contigs from these broken fragments. The purpose is to compare real scaffolds to these reconstructed contigs to gauge the usefulness of the scaffolding step. From QUAST analysis, constructed scaffolds covered 97.865% of genome. This was also the case for the broken alignment. Mismatch rate was 11.7 per 100 kilo base pairs. Additionally, there were zero misassemblies, this means that there were zero relocations, translocations and inversions. Orange blocks in the Icarus browser represented “similar misassembled blocks” [11]. This meant that they are similar because they satisfy the similarity conditions. To be labeled as similar, they must both either have correct contigs or have fragments of misassembled contigs, exceed a length threshold, and start and end within a certain range of each other. 

## Discussion

Even though QUAST reported that 97.865% of the genome was covered by this assembly, I am not comfortable calling it a finished genome. In this lab we aimed to assemble a Staphylococcus genome, a bacteria involved in antibiotic resistant infections that is used to benchmark assembly methods. Given the compact nature of bacterial genomes and their high gene density, there is a chance that in the ~2.2% of missed genome there are some genes left out critical to understanding the bacteria itself. QUAST predicted that there are 2,845 unique genes in this assembly while NCBI reports that there are actually 2,872 genes. This adds to my hesitancy to call this a completed genome because what if the gene of interest we hope to study falls in that 1% of missed genes? This concern is slightly aided by QUAST reporting that 2816 genes were captured in assembled scaffolds with only 58 genes partially captured in assembled scaffolds. 

GapCloser was able to fill 31 of the initial 94 gaps left after scaffolding. To improve gap filling and properly connect the disconnected scaffolds, I could incorporate sequence reads from alternate sequencing technologies or other sequencing runs. Given the low quality score reads, the data we started off with might not have captured these regions of the genome effectively. 

Based on the data provided by the class, as kmer size increased, both the initial max contig length assembled by minia and resulting N50 score increased as well. This makes sense because longer kmers have less of a chance to have random regions of overlap, creating more tightly formed contigs. However, there is definitely a point of diminishing returns with kmer size because it is not only limited by the sequencing technology, but also by the computational complexity of assembling longer reads. 

Raw data was much better for scaffolding. This could potentially be a result of the very low quality scores. Sickle might have over trimmed reads resulting in significantly less useful data. This left the impression that using more data is better than using the best data. Also, I expected GapClosure to result in a more significant increase in largest contig size. The maxing out of contig size could have been limited by the quality of data I began with or it might have been a technical limitation as my genome might have already been almost completely assembled. 

## Citations
1. Andrews S. [Internet]. FastQC A Quality Control tool for High Throughput Sequence Data. Babraham Bioinformatics; [cited 2018Apr10]. Available from: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ EcoliWiki [Internet].

2. John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007), DOI:10.1109/MCSE.2007.55

3. Wes McKinney. Data Structures for Statistical Computing in Python (Pandas), Proceedings of the 9th Python in Science Conference, 51-56 (2010)

4. Luo et al.: SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler. GigaScience 2012 1:18.

5. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler, QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075. doi: 10.1093/bioinformatics/btt086 First published online: February 19, 2013

6. Marten Boetzer, Christiaan V. Henkel, Hans J. Jansen, Derek Butler, Walter Pirovano; Scaffolding pre-assembled contigs using SSPACE, Bioinformatics, Volume 27, Issue 4, 15 February 2011, Pages 578–579, https://doi.org/10.1093/bioinformatics/btq683

7. Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files (Version 1.33) [Software]. Available at https://github.com/najoshi/sickle.

8. Chikhi R, Guillaume R. Space-Efficient and Exact de Bruijn Graph Representation Based on a Bloom Filter (minia). Springer, Lecture Notes In Computer Science. 2012;7534:236–48.

9. Marçais G, Kingsford C. JELLYFISH - Fast, Parallel k-mer Counting for DNA [Internet]. National Science Foundation; 2011. Available from: http://www.cbcb.umd.edu/software/jellyfish/

10. Nagarajan N, Pop M. Sequence assembly demystified. Nature Reviews Genetics. 2013;14(3):157–67.

11. Alla Mikheenko, Gleb Valin, Andrey Prjibelski, Vladislav Saveliev, Alexey Gurevich, Icarus: visualizer for de novo assembly evaluation, Bioinformatics (2016) 32 (21): 3321-3323. doi: 10.1093/bioinformatics/btw379 First published online: July 4, 2016

 

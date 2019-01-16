# Molecular Investigator: Deep Sequencing Variant Identification for Classifying Mutations in Flu Virus Allowing for Transmission to Vaccinated Hosts
##### Mauro Chavez

## Abstract
This week’s lab focused on the analysis of a “roommate’s” flu virus by deep sequencing for variants.  The roommate’s flu was able to infect someone who had previously been given a flu vaccine despite an HI test claiming the roommate’s strain was covered by the vaccine. Analysis showed it had developed quasispecies with a mutation in an epitope of the hemagglutinin gene. This meant that an altered hemaggultinin protein appeared on the virus cell membrane, making it unrecognizable by antibodies that were programmed in response to the standard H3N2 hemagglutinin protein included in this year’s vaccine.

## Introduction
Flu vaccines contain inactivated proteins from the flu virus. Once inside the body, while the inactivated viral proteins cannot cause harm, they are still recognized as foreign objects so the body proceeds to make antibodies that target regions of these proteins. Specifically, antibodies are created to target hemagglutinin proteins of the viral cells. Hemagglutinin is a viral surface protein that helps virus cells identify host cells. Epitopes are the areas of antibody-targeted molecules that are recognized and attached to by the antibody. This means that mutations in viral RNA that code for these epitopes have the potential to make resulting viral molecules unrecognizable by previously effective antibodies.  

Influenza mutations are so frequent that a single virus gives rise to many different viral species within a single host, resulting in a quasispecies.  This means that there are many rare mutations across many related viral strains increasing the possibility that one of them alters the sequence of an epitope. Existing antibodies will target viral epitopes that they recognize, creating a selective pressure towards viral cells that don’t contain these epitopes as a result of mutation. This process is called antigenic drift, where mutations in regions that code for epitopes accumulate in a viral population. 

This week involved the analysis of a flu virus sequence that was suspected to infect previously vaccinated hosts. After a hemagglutination inhibition (HI) assay, the virus seemed to match the profile of a flu strain covered by the most recent flu vaccine. However, this assay is not sensitive enough to detect rare strains in a larger viral population. In response, the flu virus’s hemagglutinin genes were deep sequenced using an Illumina single-end run and mutations were profiled to investigate if a quasispecies containing a epitope mutation existed in the population. Quasispecies represent a mixed population where mutations might occur in a very small percentage of organisms. Deep sequencing allows for many reads to be taken from a specific region of the genome, increasing the chances that these mutations are captured in the reads. Given deep sequencing reads from the various flu quasispecies, the goal was to identify the rare epitope sequence variant that resulted in hemagglutinin unrecognizable by flu vaccine created antibodies.  

This work involved the correction of next generation sequencing error that arose during the analysis procedure. During clonal amplification on Illumina sequencing chips, incorrect bases can be incorporated into sample sequence, resulting in substitution errors [7]. Because this lab involved scanning for rare variants, a data-filtering step was required to distinguish sequencing errors from real rare sequence variants. This involved sequencing the HA gene of three isogenic control samples of the standard flu genome, the reference genome mutated flu reads were mapped to. From here, any mutations detected in the control samples could be classified as sequencing errors and not true genetic variants. 
 
## Methods
In this lab, the roommate.fastq file contained deep sequencing reads of the mixed population quasispecies flu virus that was able to infect someone who had previously received the flu vaccine, specifically its hemagglutinin gene. Files SRR1705858.fastq, SRR1705859.fastq, and SRR1705860.fastq where isogenic replicates of the standard flu (H3N2 influenza) virus used for the quality control step. In sequencing, at least 151 cycles were run because the longest read length was 151 bases. However, based on the information provided in lab alone, it is not possible to tell exactly how many cycles were run. Data was pre-processed because the received reads were aligned directly to reference genome and no trimming or per base quality scores analysis were required. Additionally, across all four samples, reads were reported at various lengths, which meant low quality bases had already been trimmed. The analysis of read lengths was done using awkCommand1.txt as shown below.
```
cat roommate.fastq | awk 'NR%4==0 {print length}' | sort -n | uniq -c
```

The reference sequence for this analysis called for a standard H3N2 hemagglutinin gene.  EntrezDirect[5] was used to remotely download the sequence with id: KF848938.1 from the NCBI database using the efetch command and saved in .fasta file format.  Alignment was completed using BWA-MEM [1]. First the reference file was indexed. The alignment of roommate reads to the reference and compressing of SAM file into BAM file was completed by using bwa mem and piping the output into samtools[2] view utility and then sort utility, for all these steps default parameters were used. The resulting BAM file was then indexed using the samtools index utility. An mpileup file was created using the roommate BAM file and reference genome with –d parameter set to 1000000. This prevented samtools from ceasing piling up reads at a base after the default 8000 calls. It was important to change this parameter because otherwise, after reaching 8000 calls and stopping piling, a rare variant nucleotide could potentially be lost if it wasn’t in one of the first 8000 reads mapped to that position. 

VarScan[3] was first used to identify variants with minimum frequency 0.95 and results were outputted in .vcf file format. This meant finding variants that were common across most of the quasispecies in the roommate sample. awkCommand2.txt was then used to pull the reference base, position, and mutant base from the VarScan output.
```
cat roommate.vcf | grep -v "^#" | awk '{if (NR>24); print $2, $4, $5}'
```
Given these mutations, online sequence editor WebDSV[4] was used to view the reference genome. Mutations were recorded as which codon they altered, the resulting mutation type, and the change to the amino acid at that position in the translated protein. VarScan analysis was then repeated with a minimum variant frequency 0.001 to identify rare variants. The output of this analysis was then parsed useing awkCommand3.txt to pull mutation position, reference base, alternative base, and mutation frequency. 
```
cat roommate_rare.vcf | grep -v "^#" | awk '{print $2, $4, $5, $10}' | awk -F '[ :]' '{print $1, $2, $3, $10}'
```

In order to confidently classify rare variants as true variants and not the result of sequencing errors, the three isogenic H3N2 sequence reads had to be put through the same analysis pipeline. This process was automated using UNIX scripts. pipelineAutomation1.txt completed the aligning of each isogenic sample to the reference, BAM file creation, BAM file indexing, and read piling up. 
```
for x in SRR1705858 SRR1705859 SRR1705860
do
 echo $x
 echo "/home/linux/ieng6/cs185s/public/week2/$x.fastq"
 echo "align"
 bwa mem KF848938.1.fasta /home/linux/ieng6/cs185s/public/week2/$x.fastq | samtools view -S -b | samtools sort > $x.bam
 echo "index"
 samtools index $x.bam
 samtools mpileup -d 1000000 -f KF848938.1.fasta $x.bam > $x.mpileup
done
```
pipelineAutomation2.txt ran VarScan and parsed the resulting output files for the data of interest, primarily location, reference base, alternate base, and frequency. The exact same parameters were used here as for the analysis of the room mate sample in producing the rare variants.

```
for x in SRR1705858 SRR1705859 SRR1705860
do
 java -jar /home/linux/ieng6/cs185s/public/tools/VarScan.jar mpileup2snp $x.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > $x.vcf
 cat $x.vcf | grep -v "^#" | awk '{print $2, $4, $5, $10}' | awk -F '[ :]' '{print $1, $2, $3, $10}' > output$x.txt
done
```

True variants were then separated from sequencing errors using statistical filtering. Of the rare variants identified from the three isogenic controls, the average variant frequency and standard deviation were calculated. From there, roommate rare variants that fell within 3 standard deviations from the average control variant frequencies were removed from variant pool. The remaining ones were studied using WebDSV to record the codon change and resulting amino acid change.  These mutations were then studied under the context of the Munoz et al paper to determine if they fell within epitope regions of hemeagglutinin. Rare variant frequency data across all four samples and the calculations done for them are contained within CSE185Lab2Data.xlsx.


## Results

The results from read alignment to the reference genome are reported in tables 1A and 1B. Interestingly, when using samtools utility flagstat, the information in table 1A was reported. However, using wc –l to count the number of lines in each fastq file
and dividing by 4 provided the data in row 1 of table 1B. The number of mapped reads in table 1B was calculated using the number of initial reads minus the number of lines in each BAM file when viewed using samtools view –f4.  I am currently uncertain why this discrepancy in the number of reads exists because this method seemed to work last week. However, in both cases, the difference in initial reads and mapped reads was consistent. The difference in read count may be because reads from the fastq file were mapped more than once by samtools, increasing the effective number of total reads.

|     | roommate | SRR1705858 | SRR1705859 | SRR1705860 |
| --- | -------- | ---------- | ---------- | ---------- |
| Initial reads | 287102 | 256744 | 233451 | 250184 |
| Mapped reads | 283672 | 256658 | 233375 | 250108 |
| % mapped | 98.81% | 99.97% | 99.97% | 99.97% |

*Table 1A: sample reads, resulting number of reads mapped to reference genome, and percent of reads mapped reported by samtools flagstat*

|     | roommate | SRR1705858 | SRR1705859 | SRR1705860 |
| --- | -------- | ---------- | ---------- | ---------- |
| Initial reads | 286739 | 256586 | 233327 | 249964 |
| Mapped reads | 283309 | 256500 | 233251 | 249888 |
| % mapped | 98.80% | 99.97% | 99.97% | 99.97% |

*Table 1B: sample reads, resulting number of reads mapped to reference genome, and percent of reads mapped calculated using number of lines in fastq and output files*

Once rare variants were pulled from the output of mapping control reads to the reference genome, the average variant frequency and standard deviation for each of the three samples was calculated in excel in addition to the variant frequency and standard deviation across all three controls. This data is presented in table 2. 

|   | SRR1705858 | SRR1705859 | SRR1705860 | Total |
| - | ---------- | ---------- | ---------- | ----- |
| Avrg variant freq | 0.26% | 0.24% | 0.25% | 0.25% |
| Variant freq standard devidation | 0.072% | 0.052% | 0.078% | 0.069% |

*Table 2: data from CSE185Lab2Data.xls, average and standard deviation for variant frequency in each control sample and across each variant across all samples*

When VarScan was run for the roommate reads aligned to the reference genome with minimum variant frequency set to 0.95, the results in table 2 were returned. In table 2, the 4th through 9th columns were filled out by looking at the reference genome in WebDSV and recording the nucleotide change in each codeon and resulting amino acid change. 

| Position | REF | ALT | Original Codeon | Mutated Codeon | Original Amino Acid | Pos in Protein | Mutated Amino Acid | Mutation Type |
| -------- | --- | --- | --------------- | -------------- | ------------------- | -------------- | ------------------ | -------------- | 
| 72       |  A  |  G  | ACA | ACG | Thr | 24 | Thr | S |
| 117      |  C  |  T  | GCC | GCT | Ala | 39 | Ala | S |
| 595      |  G  |  T  | GCA | TCA | Ala | 199 | Ser | M |
| 774      |  T  |  C  | TTT | TTC | Phe | 528 | Phe | S |
| 1008     |  T  |  G  | GCT | GCG | Ala | 336 | Ala | S |
| 1260     |  A  |  C  | CTA | CTC | Leu | 420 | Leu | S |
| 1339     |  T  |  C  | TTG | CTG | Leu | 447 | Leu | S |

*Table 3: Common mutations from roommate sample identified with min variant frequency set to 0.95*

To detect rare variants, VarScan was run a second time on the piled up roommate reads with minimum variant frequency 0.001. This returned a large list of variants. From this list, variants were removed that fell within 3 standard deviations of the average variant frequency across all three isogenic controls. Variants that remained after this filtering that didn’t appear in table 3 are reported in table 4. Again, WebDSV was used to fill out the 5th through 10th columns. 

| Position | Freq | REF | ALT | Original Codeon | Mutated Codeon | Original Amino Acid | Pos in Protein | Mutated Amino Acid | Mutation Type |
| -------- | ---- | --- | --- | --------------- | -------------- | ------------------- | -------------- | ------------------ | -------------- |
|495|1.04%|C|T|AAC|AAT|Asn|165|Asn|S|
|910|0.73%|G|A|GCC|ACC|Ala|304|Asn|M|
|1293|61.82%|G|A|CTG|CTA|Leu|431|Leu|S|
|1521|1.12%|G|A|CTG|CTA|Leu|507|Leu|S|

*Table 4: Rare mutations from roomate sample identified with min variant frequency set to 0.001 that were not previously identified as common mutations*

Finally, after looking at all of the variants reported above, mutation at nucleotide 910 creating a change in the 304th amino acid falls in Epitope C of viral surface protein hemagglutinin based off findings provided by Munoz et al[6]

## Discussion
Mutations reported in both tables 3 and 4 are most likely real mutations. When rare variants were scanned from the roommate’s sequences and the filtering process was applied to them, the list of remaining variants matched the results of tables 3 and 4 combined. This means that even after filtering was applied to remove sequencing errors, these common mutations showed up both before and after filtering. Some rare mutations appeared after filtering which added to the list. 

The roommate’s flu was able to infect someone who had received the flu vaccine. The reference genome used in this analysis was a standard H3N2 influenza, one that was covered by the 2017/2018 flu vaccines. The roommate’s quasispecies flu mutation at nucleotide 910 fell in an epitope region. This means that there existed a strain of virus within the roommate that expressed a mutated hemaggultinin gene. This mutation, being in Epitope C of the HE gene as identified by Munoz et al[6], prevented previously effective antibodies from recognizing the viral cells. Since the antibodies were “programmed” using standard H3N2 influenza proteins, once the epitope in the roommates viral strain was mutated, they were no longer useful. This is how a vaccinated person was able to get the flu from the roommate, the antibodies built up as a response to the vaccine could not target the altered hemagglutinin protein. 

There are other ways of controlling error in deep sequencing experiments, Robasky et al discussed a few in their paper “The role of replicates for error mitigation in next-generation sequencing” [7]. First they mention that “filtering for sequencing read depth, base call quality, short-read alignment quality…” are “post-processing techniques [that] help to reduce uncertainty in the final genotyping variant call” [7]. These techniques were not required in this lab because the reads had been preprocessed before analysis began but they are important because simple computational steps can help highlight the important aspects of ones data. An additional way to reduce error in deep sequencing is to include “Cross-platform replicates”. Robasky et al write “each sequencing platform introduces unique biases and error types”[7]. Combining results from different procedural backgrounds helps reduce error by providing different sources for a single Variant identification. Robasky et al report that while this “greatly reduce[s] the number of false-positive variants”, it may lose many true variants in the process [7]. A third way to reduce error is by use of “Technical replicates”. Robasky defines them as “the repeat analysis of the exact same sample” [7]. Repeating the generation of data allows for the accumulation of stronger data to support one's conclusion. Being able to report with certainty that the same result can be reached through many experimental trials adds to its validity. 

Error control is a crucial step to identifying rare variants in deep sequencing data because of the nature of the variants. Rare variants may only appear in a small percentage of reads. Additionally, when dealing with virus’s that have a high replication and mutation rate, these mutations might remain in a small fraction of the quasispecies. When dealing with small frequencies of variants, the line between sequencing error and rare variant is blurred. When a rare variant has the same frequency as a sequencing error, additional steps must be taken to sperarate the two. 

There were 33 variant positions reported by VarScan that appeared in all 3 reference control sequences as identified in columns R through U in CSE185Lab2Data.xls. Alternatively to the method used in this lab, the average and standard deviation across the three reference replicates could be calculated to create a position specific filter. This data could be very useful in creating position specific error gradient where different positions in reads have different tolerance for errors. This could potentially differentiate between areas of the genome that are harder to sequence and areas that are easier to sequence. Having a different error tolerance for different positions acknowledges that as longer sequences are read, the error rate increases, and that certain areas of the genome might be less stable when put through the deep sequencing workflow. 

## Citations
1. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
2. Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
3. Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111
4. Cermak V. WebDSV [Internet]. molbiotools.com; 2014. Available from: http://www.molbiotools.com/WebDSV/index.html 
5. Kans J. Entrez Direct: E-utilities on the UNIX Command Line [Internet]. Entrez Programming Utilities Help [Internet]. U.S. National Library of Medicine; 2013 [cited 2018Apr16]. Available from: https://www.ncbi.nlm.nih.gov/books/NBK179288/ 
6. Muñoz ET, Deem MW. Epitope analysis for influenza vaccine design. Vaccine. 2005Jan19;23(9):1144–8. 
7. Robasky K, Lewis NE, Church GM. The role of replicates for error mitigation in next-generation sequencing. Nature Reviews Genetics. 2013Oct;15(1):56–62.

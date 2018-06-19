# CSE 185 Lab 6 Worksheet

##### Mauro Chavez
##### 5/14/2018

## Part 1: Mass spectrometry

**What level is the first scan?**

First scan in XML file taken at level 2

**How many peaks are in the first scan?**

205 peaks were in the first scan

**How many scans are in the entire file?**

Across the entire file, there were 21409 scans

**In the first spectrum that loaded (spectrum ID 4970), what is the m/z ratio of the most intense peak? (You can see this by hovering your mouse over that peak).** 

After loading into PRIDE inspector, spectrum 4970 had a m/z ratio of 555.065 at its most intense peak

**In the summary tab, compute the summary charts. Then: what are the two most abundant precursor ion charges?**

Precursor iron charges 2 and 3 were most abundant with 15231 and 3347 respective intensities. 

**How many spectra are in our filtered file?**

In the file filtered by msconvert, there were 1001 spectra

**The top hit is a contaminant! What is it? Give one sentence explanation of why this shows up in our list.** 

Top hit identified by mascot was a contaminant, Trypsin. However, this makes sense because Trypsin was added to samples for fragmentation

**List the top 5 scoring protein hits (not contaminants!)** 

Excluding Trypsin and other contaminants, the top 5 scoring protein hits were...

| Accession number | Score | Mass | Sig Matches |
| ---------------- | ----- | ---- | ----------- |
| NUCL_MOUSE | 161 | 76734 | 7 |
| DHX9_MOUSE | 129 | 150692 | 4 |
| TR150_MOUSE | 114 |  108171 | 5 |
| EF1A1_MOUSE | 114 | 50424 | 4 | 
| PRP8_MOUSE | 112 | 274754 | 9 | 

**How many peptides were found to match the top hit, NUCL_MOUSE (Nucleolin)?** 

7 peptides matched Nucleolin (NUCL_MOUSE)

**Briefly investigate the function of the top 3 genes. Are any known to be related to telomere function?** 

UniProt database was used to investigate protein function of top 3 hits. NUCL_MOUSE, also known as Nucleolin, is understood to initiate the decompressing of chromatin by binding to H1 histone. It may play a role in pre-rRNA transcription, transcriptional elongation, and ribosome assembly. It also binds to short RNA fragments with ` 5'-UUAGGG-3'` repeats stronger than single stranded telomeric regions do. DHX9_MOUSE is labeled as an ATP-dependent nucleic acid helicase that unwinds DNA and RNA. It is multifunctional, playing roles in DNA replication, post-transcriptional regulation of RNA, and transcriptional activation. TR150_Mouse plays a role in pre-mRNA splicing. Of the top three hits, Nucleolin was the only one associated with telomeric function.

## Part 2: De novo peptide sequencing

**How many peaks are in this file? What is the range of m/z ratios? (i.e. min and max). How does the max ion mass here compare to the precursor ion mass? What is the most abundant m/z ratio?** 

There are 1108 peaks in example_peptide.mgf with m/z ratios ranging from 279.199249 (min) to 1958.194306(max). Max ion mass is 1958.194306 which equals peptide mass (precursor ion mass). Most abundant m/z ratio was 768.640625 with intensity 7950.795898.

**Include your table in your worksheet. Explain your answer for at least two more amino acids.** 

| Amino acid | Position (b-series) | Position (y-series) | Peak (b) | Peak (y) |
|--------|---------|---------|-------|------|
| I | 1 | 18 | 114.0919 | 1845.923 |
| A | 2 | 17 | 185.0371 | 1774.89|
| ? | 3 | 16 | | |
| ? | 4 | 15 | | |
| ? | 5 | 14 | | |
| ? | 6 | 13 | | |
| ? | 7 | 12 | | |
| ? | 8 | 11 | | |
| ? | 9 | 10 | | |
| ? | 10 | 9 | | |
| ? | 11 | 8 | | |
| ? | 12 | 7 | | |
| ? | 13 | 6 | | |
| ? | 14 | 5 | | |
| ? | 15 | 4 | 1507.841 | 450.226084 |
| T | 16 | 3 | 1654.91092 | 303.15767 |
| G | 17 | 2 | 1755.9586 | 202.119 |
| K | 18 | 1 | 1812.90 | 145.0977 |

- To identify amino acid A, I noticed range of high intensity peaks in 1774.0265 to 1776.111816 region. The distance from 1812 to 1775 (37) is too small to represent an amino acid so it must correspond to one at the front. The difference between 1845 and 1774 is ~71.9 which is very close to 71.037, the mass of Serine. 
- To identify amino acid G, I noticed a range of peaks in 1752-1758 region. 1812-1756 = 56.98 which represents the mass shift provided by a Glycerine. 
- There were large peaks in 1653 to 1658 range, subtracting this range from 1755 gave masses close to Theronine's 101.  
- I saw large peaks in 451.382812 to 456.441742 region. The largest peak being at 452.38281. 303.151767 - 452.38281 = 149.23, which closely matches the mass of Phenylalanine.

**This was pretty tedious! If you were going to write a program to automate this process, how would you do it? (max 2-3 sentences about your idea.** 

To automate this process I would make a graph where nodes are peaks and edges connect peaks if the difference between them is close to an amino acid weight. Finding a path from the start node with mass zero to the node that represents the mass of the peptide would represent the string of amino acids in the peptide. 

**What is the top scoring sequence? How close is it to your guess? Note: neither answer is probably completely correct!**

The top scoring sequence from pepnovo was `PQELLSQLVQYTGK` which matches the last three amino acids I guessed but doesnt match the first two. 

**Hypothesize why some peaks are missing from our spectrum file.** 

Based off the theoretical spectrum for `PQELLSQLVQYTGK`, smaller mass peaks are missing from our spectrum, possibly because the device used fore measurements was not sensitive enough. Some peaks do appear in the spectrum (1717 in range 17.16 to 1719) but the intensity is not large enough to make them significant next to some of the neighboring peaks. This is a result of the noisy and messy spectrum we were provided. 786.64 is a very large peak that shows up in spectrum but not in the theoretical one, its inclusion might be skewing all calculations made after it is incorporated. Alternatively, it might be a signal given off by Trypsin.  

**Suggest a way we could modify our mass spec experiment to distinguish between amino acids with near identical masses, such as Isoleucine and Leucine.** 

To solve the issue of Isolucine and Leucine masses being identical, the molecule of interest could be preprocessed changing the structure of one of these two amino only by adding or removing a component. Alternatively, an enzyme could be created that binds to Isolucine or Lucine adding weight, giving it a more distinct mass footprint.

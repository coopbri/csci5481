### CSCI 5481 Homework 4
###### Brian Cooper

## Code Explanation
I used my code (in `homework4.py`) for all four of the steps in the homework.

For the phylogeny construction and visualization in step 4, I first used a sequence tool conversion tool at http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_nexus.php to convert the FASTA files to <a href="https://en.wikipedia.org/wiki/Nexus_file">NEXUS</a>-format files. Then, I used the phylogeny visualization tool at http://www.phylogeny.fr/simple_phylogeny.cgi to generate visual trees based on the NEXUS files.

## Guide to Deliverables
- _Source files (any code that you used for Step 1, 2, 3, 4)_
  - `homework4.py`
- _Readme file explaining how you used your code (text)_
  - `README.md` (this file)
- _Step 1: File giving start and end position of each variable region_
  - `regions.txt`
    - Also, a visual: `regions.png`
- _Step 2: File containing identity values_
  - `variability.txt`
- _Step 3: Figure showing plot of identity values._
  - `variability.png`
- _Step 4: Fasta files containing variable regions 1 and 4. Figures of the V1 tree, V4 tree, and whole-16S tree. Brief text response to questions._
  - FASTA region files: `r1.fna`, `r4.fna`, `whole.nexus`
    - Intermediate files: `r1.nexus`, `r4.nexus`, `whole.nexus`
  - Figure files: `v1.png`, `v4.png`, `whole.png`
  - Text responses: `responses.md`

<center>
  <h1>Reading Response | Metagenomics - DADA2</h1>
</center>

### Key points (1-3 sentences)
A software called DADA2, which extends the original DADA algorithm by providing more robust species and taxonomic assignment/inference, is introduced in this paper. With its introduction, DADA2 was tested on various communities of data, including a "Balanced" set (bacteria and archaea from various habitats), an "HMP" set (well-separated strains, mainly from the human body), and an "Extreme" set (bacteria from human gastrointestinal tract). With these test sets, DADA2 was shown to perform quite well compared to various other software (such as UPARSE and MED) -- DADA2 made many exact output sequences on the test data sets.

### Added to reading (1-2 sentences)
The authors mention runtime analysis, but don't appear to go into detail about the analysis itself (DADA2 and how it compares to other packages). A more clear position of DADA2's time and space complexity would augment the reading well.

### Question(s)
The nucleotide variance difference detection is impressive. What other algorithms perform well for difference resolution?

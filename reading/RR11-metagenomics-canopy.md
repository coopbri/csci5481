<center>
  <h1>Reading Response | Metagenomics - Canopy</h1>
</center>

### Key points (1-3 sentences)
In this paper, the authors suggest a metagenomics processing method for microbial genomic assembly and identification of microbial entities that involves using sets of dense samples of similar types (and further avoiding the requirement of reference sequences). To do this without reference sequences, they propose selecting a "seed" gene (randomly), and place that gene with other genes of similar abundance profiles, forming what they call a "canopy" (which specifically involves genes with a Pearson correlation coefficient of 0.9 or greater). Metagenomic species (MGS) and co-abundance gene groups (CAGs) were used as groupings for the MGS canopy process.

### Added to reading (1-2 sentences)
I tried accessing the source code from their Git repo link (http://git.dworzynski.eu/mgs-canopy-algorithm), but the link leads to Bitbucket (must be their old online repo platform) and the bare domain leads to a DreamHost landing page (so they are paying for hosting but the website appears to be unused). I did find the source code on GitHub so I assume they moved it to that, but the reading would benefit from an updated link (despite the paper being already printed, of course).

### Question(s)
The authors mention a method of reducing "chimeric assemblies of closely related strains" by using only reads from a single sample for each chimeric assembly in question. This seems fine to me, but what do the authors mean when they follow this information with the statement, ">100 MGS could be assembled from multiple samples, hence, the total number of high-quality draft genome assemblies was 360"?

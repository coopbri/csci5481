<center>
  <h1>Reading Response | Metagenomics - Kraken</h1>
</center>

### Key points (1-3 sentences)
This paper introduces the now-popular taxonomic sequence classification tool called Kraken. Kraken is most often used in metagenomics sequence labeling tasks, and is demonstrated to be exceptionally fast and accurate, largely due in part for its exact-database-matching of k-mers. The core drawback of Kraken is that it requires large amounts of memory (70GB+ for the "default database" setup), although the authors provide a method of k-mer removal that drastically lowers the memory requirement.

### Added to reading (1-2 sentences)
I really like how this paper set up the breadth of the metagenomics field and explained why Kraken is useful before proposing the tool itself in the Background section. I think the reading could benefit from deeper elaboration of their kernel space computer memory allocation and management, which leads into my question below.


### Question(s)
The authors mentioned that they normalized the tests by reading the database files into memory three times to have this content in the OS cache. I am curious about how they measured this (at the OS/kernel level) and why they chose three times.

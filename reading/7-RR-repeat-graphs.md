<center>
  <h1>Reading Response | Assembly of ... Repeat Graphs</h1>
</center>

### Key points (1-3 sentences)
This article introduced the Flye algorithm for genome assembly and compared it to several "state-of-the art" assemblers: Canu, Falcon, HINGE, Miniasm, and MaSuRCA. Flye generates a repeat graph to store information about repeating segments in a genome for unbridged repeats (repeats are "bridged" if 1+ read(s) contains one of its copies). Flye was shown to assemble a larger amount of the human genome than Canu and MaSuRCA, assemble more accurately than HINGE on complex eukaryotes (YEAST and WORM were tested),  was faster than Falcon on assembling the WORM dataset, and had a higher sequence identity than Miniasm.

### Added to reading (1-2 sentences)
The reading appears to hastily end with saying that there is potential room for growth in research for assembly algorithms that handle unbridged repeats. I think some more direction in this regard would augment the reading well (beyond their explanation of low number of segmental duplications).

### Question(s)
I think this is from a lack of understanding about these two algorithms, but the paper mentioned that "HINGE does not distinguish between complete and semi-complete assemblies". I could not find this in the HINGE introduction paper (https://genome.cshlp.org/content/27/5/747.full) and was wondering how the authors determined this.

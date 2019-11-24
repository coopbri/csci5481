<center>
  <h1>Reading Response | RNA-Seq TopHat</h1>
</center>

### Key points (1-3 sentences)
This paper covers the (at-the-time) limitations of identifying novel splice junctions with RNA-Seq. Following this, TopHat is introduced, which is an algorithm for RNA-Seq read alignment without reliance on such (known) splice junctions. This quality makes TopHat useful (not necessarily efficient, but certainly useful).

### Added to reading (1-2 sentences)
The statistic mentioned in the "Methods" section (Dij) was confusing to me. I think more explanation of the statistic formula would be useful.

### Question(s)
The paper mentions "coding regions must still be distinguished from UTRs and non-coding RNAs." In what ways has this task been improved since this paper's publishing? Similarly, TopHat appears to be fast compared to other approaches at the time (with 2.2M reads per CPU hour), but how does it hold up with modern techniques?

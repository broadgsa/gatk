## Notes on downsampling in HC/M2

http://gatkforums.broadinstitute.org/gatk/discussion/8028/notes-on-downsampling-in-hc-m2

<p><strong>This document aims to record some developer notes for posterity. Contents were generated July 24, 2015 and are not guaranteed to be up to date. No support guarantee either.</strong></p>
<hr />
<h3>Arguments and Parameters</h3>
<ul>
<li>&quot;@Downsample&quot; annotation in the class definition for HC/M2 controls the height of the pileup used for active region determination -- HC has default coverage 500, M2 downsampling takes on ActiveRegionWalker default, which is 1000</li>
<li>&quot;@ActiveRegionTraversalParameters&quot; argument controls the maximum number of reads that can possibly be processed; has default maxReadsToHoldInMemoryPerSample() = 30,000 and across all samples maxReadsToHoldTotal() = 10,000,000 -- these are not currently overridden in HC or M2</li>
<li>maxReadsInRegionPerSample and minReadsPerAlignmentStart (arguments in HC, hard-coded in M2 right now) loosely control the number of reads that go into the assembly step; default is 10K and 10 for HC, hard-coded 1000 and 5 for M2</li>
</ul>
<h3>Relevant Code</h3>
<ul>
<li>TraverseActiveRegions.java does a lot of the data management (including some downsampling) and does the iteration over ActiveRegions in traverse()</li>
<li>TraverseActiveRegions takes in the ActiveRegionTraversalParameters annotations and creates a TAROrderedReadCache (this is where all the reads that get passed to the Walker are stored)</li>
<li>TAROrderedReadCache contains a ReservoirDownsampler </li>
<li>ReservoirDownsampler is unbiased with regard to read start position; gets initialized with maxCapacity = min(maxReadsToHoldTotal, maxReadsToHoldInMemoryPerSample*nSamples)</li>
<li>Reads that go into the Walker's map() function get downsampled by the ReservoirDownsampler to exactly maxCapacity if they exceed the maxCapacity -- at this point this is the most reads you can ever use for calculations</li>
<li>Reads that go into the assembly step (already filtered for MQ) get downsampled by the LevelingDownsampler to approximately maxReadsInRegionPerSample if the number of reads exceeds maxReadsInRegionPerSample
<ul>
<li>(my maxReadsInRegionPerSample is one step-through was 1037, but my downsampleReads was 3003 over 100bp, so it seems pretty approximate)</li>
</ul></li>
<li>The LevelingDownsampler is intentionally biased because it maintains a minimum coverage at each base as specified by minReadsPerAlignmentStart</li>
<li>ActiveRegionTrimmer.Result trimmingResult in the Walker's map() function recovers reads (up to theTAROrderedReadCache maxCapacity) by pulling them from the originalActiveRegion, but trims them to variation events found in the (potentially downsampled) assembly</li>
<li>Genotyping is performed based largely on the set of reads going into the map() function (M2 filters for quality with filterNonPassingReads before genotyping)</li>
</ul>
<h3>Worst Case M2 Behavior</h3>
<ul>
<li>Highest coverage M2 call on CRSP NA12878 SM-612V4.bam vs SM-612V3.bam normal-normal pair occurs at 7:100645781 with 4000-5000X coverage, also coverage can exceed 7000X in other BAMs</li>
<li>A lot of exons have coverage exceeding the 1000X cutoff for ActiveRegion determination with isActive(), but even with downsampling to 1000X we should still trigger ActiveRegions down to around allele fraction of ~0.8% for 4000X</li>
<li>Even the highest coverage exon in the CRSP NA12878 normal-normal calling doesn't exceed the default limit for the ReservoirDownsampler (i.e. all reads will have the potential to get genotyped)</li>
<li>In this super high coverage exon, reads are getting downsampled to ~3000 before they go into the assembly (again, controlled by maxReadsInRegionPerSample and minReadsPerAlignmentStart)
<ul>
<li>Here that's a retention of about 12.5% of reads, which seems pretty aggressive</li>
<li>The maxReadsInRegionPerSample value is 10% of what it is for HC</li>
<li>Increasing maxReadsInRegionPerSample for M2 may increase sensitivity (although honestly not based on my LUAD comparison vs. M1) but will drastically increase assembly time</li>
</ul></li>
<li>All reads that pass quality filters are genotyped according to the variants found using the downsampled assembly set</li>
</ul>
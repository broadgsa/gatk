## Downsampling

http://gatkforums.broadinstitute.org/gatk/discussion/1323/downsampling

<h4>Downsampling is a process by which read depth is reduced, either at a particular position or within a region.</h4>
<p>Normal sequencing and alignment protocols can often yield pileups with vast numbers of reads aligned to a single section of the genome in otherwise well-behaved datasets.  Because of the frequency of these 'speed bumps', the GATK now downsamples pileup data unless explicitly overridden.</p>
<p>Note that there is also a proportional &quot;downsample to fraction&quot; mechanism that is mostly intended for testing the effect of different overall coverage means on analysis results.</p>
<p>See below for details of how this is implemented and controlled in GATK.</p>
<hr />
<h2>1. Downsampling to a target coverage</h2>
<p>The principle of this downsampling type is to downsample reads to a given capping threshold coverage. Its purpose is to get rid of excessive coverage, because above a certain depth, having additional data is not informative and imposes unreasonable computational costs. The downsampling process takes two different forms depending on the type of analysis it is used with. For locus-based traversals (LocusWalkers like UnifiedGenotyper and ActiveRegionWalkers like HaplotypeCaller), downsample_to_coverage controls the maximum depth of coverage at each locus. For read-based traversals (ReadWalkers like BaseRecalibrator), it controls the maximum number of reads sharing the same alignment start position. For ReadWalkers you will typically need to use much lower dcov values than you would with LocusWalkers to see an effect. Note that this downsampling option does not produce an unbiased random sampling from all available reads at each locus: instead, the primary goal of the to-coverage downsampler is to maintain an even representation of reads from all alignment start positions when removing excess coverage. For a truly unbiased random sampling of reads, use -dfrac instead. Also note that the coverage target is an approximate goal that is not guaranteed to be met exactly: the downsampling algorithm will under some circumstances retain slightly more or less coverage than requested.</p>
<h3>Defaults</h3>
<p>The GATK's default downsampler (invoked by <code>-dcov</code>) exhibits the following properties:</p>
<ul>
<li>The downsampler treats data from each sample independently, so that high coverage in one sample won't negatively impact calling in other samples.  </li>
<li>The downsampler attempts to downsample uniformly across the range spanned by the reads in the pileup.  </li>
<li>The downsampler's memory consumption is proportional to the sampled coverage depth rather than the full coverage depth.</li>
</ul>
<p>By default, the downsampler is limited to 1000 reads per sample.  This value can be adjusted either per-walker or per-run.</p>
<h3>Customizing</h3>
<p>From the command line:</p>
<ul>
<li>To disable the downsampler, specify <code>-dt NONE</code>.  </li>
<li>To change the default coverage per-sample, specify the desired coverage to the <code>-dcov</code> option.</li>
</ul>
<p>To modify the walker's default behavior:</p>
<ul>
<li>Add the @Downsample interface to the top of your walker.  Override the downsampling type by changing the <code>by=&lt;value&gt;</code>.  Override the downsampling depth by changing the <code>toCoverage=&lt;value&gt;</code>.</li>
</ul>
<h3>Algorithm details</h3>
<p>The downsampler algorithm is designed to maintain uniform coverage while preserving a low memory footprint in regions of especially deep data. Given an already established pileup, a single-base locus, and a pile of reads with an alignment start of single-base locus + 1, the outline of the algorithm is as follows:</p>
<p>For each sample:</p>
<ul>
<li>Select <sample size> reads with the next alignment start.  </li>
<li>While the number of existing reads + the number of incoming reads is greater than the target sample size:</li>
</ul>
<p>Now walk backward through each set of reads having the same alignment start.  If the count of reads having the same alignment start is &gt; 1, throw out one randomly selected read.</p>
<ul>
<li>If we have n slots available where n is &gt;= 1, randomly select n of the incoming reads and add them to the pileup.  </li>
<li>Otherwise, we have zero slots available.  Choose the read from the existing pileup with the least alignment start.  Throw it out and add one randomly selected read from the new pileup.</li>
</ul>
<hr />
<h2>2. Downsampling to a fraction of the coverage</h2>
<p>Reads will be downsampled so the specified fraction remains; e.g. if you specify -dfrac 0.25, three-quarters of the reads will be removed, and the remaining one quarter will be used in the analysis. This method of downsampling is truly unbiased and random. It is typically used to simulate the effect of generating different amounts of sequence data for a given sample. For example, you can use this in a pilot experiment to evaluate how much target coverage you need to aim for in order to obtain enough coverage in all loci of interest.</p>
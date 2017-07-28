## Accessing reads: AlignmentContext and ReadBackedPileup

http://gatkforums.broadinstitute.org/gatk/discussion/1322/accessing-reads-alignmentcontext-and-readbackedpileup

<h3>1. Introduction</h3>
<p>The AlignmentContext and ReadBackedPileup work together to provide the read data associated with a given locus.  This section details the tools the GATK provides for working with collections of aligned reads.</p>
<h3>2. What are read backed pileups?</h3>
<p>Read backed pileups are objects that contain all of the reads and their offsets that &quot;pile up&quot; at a locus on the genome.  They are the basic input data for the GATK LocusWalkers, and underlie most of the locus-based analysis tools like the recalibrator and SNP caller.  Unfortunately, there are many ways to view this data, and version one grew unwieldy trying to support all of these approaches.   Version two of the ReadBackedPileup presents a consistent and clean interface for working pileup data, as well as supporting the <code>iterable()</code> interface to enable the convenient <code>for ( PileupElement p : pileup )</code> for-each loop support.</p>
<h3>3. How do I get a ReadBackedPileup and/or how do I create one?</h3>
<p>The best way is simply to grab the pileup (the underlying representation of the locus data) from your <code>AlignmentContext</code> object in <code>map</code>:</p>
<pre><code class="pre_md">public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    ReadBackedPileup pileup = context.getPileup();</code class="pre_md"></pre>
<p>This aligns your calculations with the GATK core infrastructure, and avoids any unnecessary data copying from the engine to your walker.</p>
<h4>If you are trying to create your own, the best constructor is:</h4>
<pre><code class="pre_md">public ReadBackedPileup(GenomeLoc loc, ArrayList&lt;PileupElement&gt; pileup )</code class="pre_md"></pre>
<p>requiring only a list, in order of read / offset in the pileup, of PileupElements.</p>
<h4>From List<SAMRecord> and List<Offset></h4>
<p>If you happen to have lists of SAMRecords and integer offsets into them you can construct a <code>ReadBackedPileup</code> this way:</p>
<pre><code class="pre_md">public ReadBackedPileup(GenomeLoc loc, List&lt;SAMRecord&gt; reads, List&lt;Integer&gt; offsets )</code class="pre_md"></pre>
<h3>4. What's the best way to use them?</h3>
<h4>Best way if you just need reads, bases and quals</h4>
<pre><code class="pre_md">for ( PileupElement p : pileup ) {
  System.out.printf("%c %c %d%n", p.getBase(), p.getSecondBase(), p.getQual());
  // you can get the read itself too using p.getRead()
}</code class="pre_md"></pre>
<p>This is the most efficient way to get data, and should be used whenever possible.</p>
<h4>I just want a vector of bases and quals</h4>
<p>You can use:</p>
<pre><code class="pre_md">public byte[] getBases()
public byte[] getSecondaryBases()
public byte[] getQuals()</code class="pre_md"></pre>
<p>To get the bases and quals as a <code>byte[]</code> array, which is the underlying base representation in the SAM-JDK.</p>
<h4>All I care about are counts of bases</h4>
<p>Use the follow function to get counts of A, C, G, T in order: </p>
<pre><code class="pre_md">public int[] getBaseCounts()</code class="pre_md"></pre>
<p>Which returns a <code>int[4]</code> vector with counts according to <code>BaseUtils.simpleBaseToBaseIndex</code> for each base.</p>
<h4>Can I view just the reads for a given sample, read group, or any other arbitrary filter?</h4>
<p>The GATK can very efficiently stratify pileups by sample, and less efficiently stratify by read group, strand, mapping quality, base quality, or any arbitrary filter function.  The sample-specific functions can be called as follows:</p>
<pre><code class="pre_md">pileup.getSamples();
pileup.getPileupForSample(String sampleName);</code class="pre_md"></pre>
<p>In addition to the rich set of filtering primitives built into the <code>ReadBackedPileup</code>, you can supply your own primitives by implmenting a PileupElementFilter:</p>
<pre><code class="pre_md">public interface PileupElementFilter {
    public boolean allow(final PileupElement pileupElement);
}</code class="pre_md"></pre>
<p>and passing it to <code>ReadBackedPileup</code>'s generic filter function:</p>
<pre><code class="pre_md">public ReadBackedPileup getFilteredPileup(PileupElementFilter filter);</code class="pre_md"></pre>
<p>See the <code>ReadBackedPileup</code>'s java documentation for a complete list of built-in filtering primitives.</p>
<h4>Historical: StratifiedAlignmentContext</h4>
<p>While <code>ReadBackedPileup</code> is the preferred mechanism for aligned reads, some walkers still use the <code>StratifiedAlignmentContext</code> to carve up selections of reads.  If you find functions that you require in <code>StratifiedAlignmentContext</code> that seem to have no analog in <code>ReadBackedPileup</code>, please let us know and we'll port the required functions for you.</p>
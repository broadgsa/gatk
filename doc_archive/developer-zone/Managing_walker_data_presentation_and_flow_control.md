## Managing walker data presentation and flow control

http://gatkforums.broadinstitute.org/gatk/discussion/1351/managing-walker-data-presentation-and-flow-control

<p>The primary goal of the GATK is to provide a suite of small data access patterns that can easily be parallelized and otherwise externally managed.  As such, rather than asking walker authors how to iterate over a data stream, the GATK asks the user how data should be presented.  </p>
<h2>Locus walkers</h2>
<p>Walk over the data set one location (single-base locus) at a time, presenting all overlapping reads, reference bases, and reference-ordered data.</p>
<h3>1. Switching between covered and uncovered loci</h3>
<p>The <code>@By</code> attribute can be used to control whether locus walkers see all loci or just covered loci.  To switch between viewing all loci and covered loci, apply one of the following attributes:</p>
<pre><code class="pre_md">@By(DataSource.REFERENCE)
@By(DataSource.READS)</code class="pre_md"></pre>
<h3>2. Filtering defaults</h3>
<p>By default, the following filters are automatically added to every locus walker.</p>
<ul>
<li>Reads with nonsensical alignments</li>
<li>Unmapped reads</li>
<li>Non-primary alignments.</li>
<li>Duplicate reads.</li>
<li>Reads failing vendor quality checks.</li>
</ul>
<h2>ROD walkers</h2>
<p>These walkers walk over the data set one location at a time, but only those locations covered by reference-ordered data.  They are essentially a special case of locus walkers. ROD walkers are read-free traversals that include operate over Reference Ordered Data and the reference genome <strong>at sites where there is ROD information</strong>.  They are geared for high-performance traversal of many RODs and the reference such as VariantEval and CallSetConcordance.  Programmatically they are nearly identical to <code>RefWalkers&lt;M,T&gt;</code> traversals with the following few quirks.</p>
<h3>1. Differences from a RefWalker</h3>
<ul>
<li>
<p>RODWalkers are only called at sites where there is at least one non-interval ROD bound.  For example, if you are exploring dbSNP and some GELI call set, the map function of a RODWalker will be invoked at all sites where there is a dbSNP record or a GELI record.</p>
</li>
<li>
<p>Because of this skipping RODWalkers receive a context object where the number of reference skipped bases between map calls is provided: </p>
<p>nSites += context.getSkippedBases() + 1; // the skipped bases plus the current location</p>
</li>
</ul>
<p>In order to get the final count of skipped bases at the end of an interval (or chromosome) the map function is called one last time with null <code>ReferenceContext</code> and <code>RefMetaDataTracker</code> objects.  The alignment context can be accessed to get the bases skipped between the last (and final) ROD and the end of the current interval. </p>
<h3>2. Filtering defaults</h3>
<p>ROD walkers inherit the same filters as locus walkers:</p>
<ul>
<li>Reads with nonsensical alignments</li>
<li>Unmapped reads</li>
<li>Non-primary alignments.</li>
<li>Duplicate reads.</li>
<li>Reads failing vendor quality checks.</li>
</ul>
<h3>3. Example change over of VariantEval</h3>
<p>Changing to a RODWalker is very easy -- here's the new top of VariantEval, changing the system to a <code>RodWalker</code> from its old <code>RefWalker</code> state:</p>
<pre><code class="pre_md">//public class VariantEvalWalker extends RefWalker&lt;Integer, Integer&gt; {
public class VariantEvalWalker extends RodWalker&lt;Integer, Integer&gt; {</code class="pre_md"></pre>
<p>The map function must now capture the number of skipped bases and protect itself from the final interval map calls:</p>
<pre><code class="pre_md">public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
    nMappedSites += context.getSkippedBases();

    if ( ref == null ) { // we are seeing the last site
        return 0;
    }

    nMappedSites++;</code class="pre_md"></pre>
<p>That's all there is to it!</p>
<h3>4. Performance improvements</h3>
<p>A ROD walker can be very efficient compared to a RefWalker in the situation where you have sparse RODs. Here is a comparison of ROD vs. Ref walker implementation of VariantEval:</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;">RODWalker</th>
<th style="text-align: left;">RefWalker</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">dbSNP and 1KG Pilot 2 SNP calls on chr1</td>
<td style="text-align: left;">164u (s)</td>
<td style="text-align: left;">768u (s)</td>
</tr>
<tr>
<td style="text-align: left;">Just 1KG Pilot 2 SNP calls on chr1</td>
<td style="text-align: left;">54u (s)</td>
<td style="text-align: left;">666u (s)</td>
</tr>
</tbody>
</table>
<h2>Read walkers</h2>
<p>Read walkers walk over the data set one read at a time, presenting all overlapping reference bases and reference-ordered data.</p>
<h3>Filtering defaults</h3>
<p>By default, the following filters are automatically added to every read walker.</p>
<ul>
<li>Reads with nonsensical alignments</li>
</ul>
<h2>Read pair walkers</h2>
<p>Read pair walkers walk over a queryname-sorted BAM, presenting each mate and its pair.  No reference bases or reference-ordered data are presented.</p>
<h3>Filtering defaults</h3>
<p>By default, the following filters are automatically added to every read pair walker.</p>
<ul>
<li>Reads with nonsensical alignments</li>
</ul>
<h2>Duplicate walkers</h2>
<p>Duplicate walkers walk over a read and all its marked duplicates.  No reference bases or reference-ordered data are presented.</p>
<h3>Filtering defaults</h3>
<p>By default, the following filters are automatically added to every duplicate walker.</p>
<ul>
<li>Reads with nonsensical alignments</li>
<li>Unmapped reads</li>
<li>Non-primary alignments</li>
</ul>
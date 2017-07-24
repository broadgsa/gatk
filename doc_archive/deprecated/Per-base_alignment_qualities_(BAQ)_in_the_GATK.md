## Per-base alignment qualities (BAQ) in the GATK

http://gatkforums.broadinstitute.org/gatk/discussion/1326/per-base-alignment-qualities-baq-in-the-gatk

<h3>This article is out of date and no longer applicable. BAQs are no longer used in GATK.</h3>
<hr />
<h3>1. Introduction</h3>
<p>The GATK provides an implementation of the Per-Base Alignment Qualities (BAQ) developed by Heng Li in late 2010.  See <a href="http://samtools.sourceforge.net/mpileup.shtml">this SamTools page</a> for more details.</p>
<hr />
<h3>2. Using BAQ</h3>
<p>The BAQ algorithm is applied by the GATK engine itself, which means that all GATK walkers can potentially benefit from it.  By default, BAQ is OFF, meaning that the engine will not use BAQ quality scores at all. </p>
<p>The GATK engine accepts the argument <code>-baq</code> with the following <code>enum</code> values: </p>
<pre><code class="pre_md">public enum CalculationMode {
    OFF,                        // don't apply a BAQ at all, the default
    CALCULATE_AS_NECESSARY,     // do HMM BAQ calculation on the fly, as necessary, if there's no tag
    RECALCULATE                 // do HMM BAQ calculation on the fly, regardless of whether there's a tag present
}</code class="pre_md"></pre>
<p>If you want to enable BAQ, the usual thing to do is <code>CALCULATE_AS_NECESSARY</code>, which will calculate BAQ values if they are not in the <code>BQ</code> read tag.  If your reads are already tagged with <code>BQ</code> values, then the GATK will use those.  <code>RECALCULATE</code> will always recalculate the <code>BAQ</code>, regardless of the tag, which is useful if you are experimenting with the gap open penalty (see below).</p>
<p>If you are really an expert, the GATK allows you to specify the BAQ gap open penalty (<code>-baqGOP</code>) to use in the HMM.  This value should be 40 by default, a good value for whole genomes and exomes for highly sensitive calls. However, if you are analyzing exome data only, you may want to use 30, which seems to result in more specific call set.  We continue to play with these values some.  Some walkers, where BAQ would corrupt their analyses, forbid the use of BAQ and will throw an exception if <code>-baq</code> is provided.</p>
<hr />
<h3>3. Some example uses of the BAQ in the GATK</h3>
<ul>
<li>
<p>For UnifiedGenotyper to get more specific SNP calls.</p>
</li>
<li>
<p>For PrintReads to write out a BAM file with BAQ tagged reads</p>
</li>
<li>For TableRecalibrator or IndelRealigner to write out a BAM file with BAQ tagged reads.  Make sure you use <code>-baq RECALCULATE</code> so the engine knows to recalculate the BAQ after these tools have updated the base quality scores or the read alignments.  Note that both of these tools will not use the BAQ values on input, but will write out the tags for analysis tools that will use them.</li>
</ul>
<p>Note that some tools should not have BAQ applied to them.</p>
<p>This last option will be a particularly useful for people who are already doing base quality score recalibration.  Suppose I have a pipeline that does:</p>
<pre><code class="pre_md">RealignerTargetCreator
IndelRealigner

BaseRecalibrator
PrintReads (with --BQSR input)

UnifiedGenotyper</code class="pre_md"></pre>
<p>A highly efficient BAQ extended pipeline would look like</p>
<pre><code class="pre_md">RealignerTargetCreator
IndelRealigner // don't bother with BAQ here, since we will calculate it in table recalibrator

BaseRecalibrator
PrintReads (with --BQSR input) -baq RECALCULATE // now the reads will have a BAQ tag added.  Slows the tool down some

UnifiedGenotyper -baq CALCULATE_AS_NECESSARY // UG will use the tags from TableRecalibrate, keeping UG fast</code class="pre_md"></pre>
<hr />
<h3>4. BAQ and walker control</h3>
<p>Walkers can control via the <code>@BAQMode</code> annotation how the BAQ calculation is applied.  Can either be as a tag, by overwriting the qualities scores, or by only returning the baq-capped qualities scores.  Additionally, walkers can be set up to have the BAQ applied to the incoming reads (<code>ON_INPUT</code>, the default), to output reads (<code>ON_OUTPUT</code>), or <code>HANDLED_BY_WALKER</code>, which means that calling into the BAQ system is the responsibility of the individual walker.</p>
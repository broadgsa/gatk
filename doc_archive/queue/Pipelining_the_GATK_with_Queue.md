## Pipelining the GATK with Queue

http://gatkforums.broadinstitute.org/gatk/discussion/1310/pipelining-the-gatk-with-queue

<h3>1. Introduction</h3>
<p>As mentioned in the introductory materials, the core concept behind the GATK tools is the walker. The Queue scripting framework contains several mechanisms which make it easy to chain together GATK walkers.</p>
<h3>2. Authoring walkers</h3>
<p>As part of authoring your walker there are several Queue behaviors that you can specify for [QScript]() authors using your particular walker.</p>
<h4>Specifying how to partition</h4>
<p>Queue can significantly speed up generating walker outputs by passing different instances of the GATK the same BAM or VCF data but specifying different regions of the data to analyze. After the different instances output their individual results Queue will gather the results back to the original output path requested by QScript.</p>
<p>Queue limits the level it will split genomic data by examining the <code>@PartitionBy()</code> annotation for your walker which specifies a <code>PartitionType</code>. This table lists the different partition types along with the default partition level for each of the different walker types.</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">PartitionType</th>
<th style="text-align: left;">Default for Walker Type</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Example Intervals</th>
<th style="text-align: left;">Example Splits</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">PartitionType.CONTIG</td>
<td style="text-align: left;">Read walkers</td>
<td style="text-align: left;">Data is grouped together so that all genomic data from the same contig is never presented to two different instances of the GATK.</td>
<td style="text-align: left;">original: chr1:10-11, chr2:10-20, chr2:30-40, chr2:50-60, chr3:10-11</td>
<td style="text-align: left;">split 1: chr1:10-11, chr2:10-20, chr2:30-40, chr2:50-60; split 2:chr3:10-11</td>
</tr>
<tr>
<td style="text-align: left;">PartitionType.INTERVAL</td>
<td style="text-align: left;">(none)</td>
<td style="text-align: left;">Data is split down to the interval level but never divides up an explicitly specified interval. If no explicit intervals are specified in the QScript for the GATK then this is effectively the same as splitting by contig.</td>
<td style="text-align: left;">original: chr1:10-11, chr2:10-20, chr2:30-40, chr2:50-60, chr3:10-11</td>
<td style="text-align: left;">split 1: chr1:10-11, chr2:10-20, chr2:30-40; split 2: chr2:50-60, chr3:10-11</td>
</tr>
<tr>
<td style="text-align: left;">PartitionType.LOCUS</td>
<td style="text-align: left;">Locus walkers, ROD walkers</td>
<td style="text-align: left;">Data is split down to the locus level possibly dividing up intervals.</td>
<td style="text-align: left;">original: chr1:10-11, chr2:10-20, chr2:30-40, chr2:50-60, chr3:10-11</td>
<td style="text-align: left;">split 1: chr1:10-11, chr2:10-20, chr2:30-35; split 2: chr2:36-40, chr2:50-60, chr3:10-11</td>
</tr>
<tr>
<td style="text-align: left;">PartitionType.NONE</td>
<td style="text-align: left;">Read pair walkers, Duplicate walkers</td>
<td style="text-align: left;">The data cannot be split and Queue must run the single instance of the GATK as specified in the QScript.</td>
<td style="text-align: left;">original: chr1:10-11, chr2:10-20, chr2:30-40, chr2:50-60, chr3:10-11</td>
<td style="text-align: left;">no split: chr1:10-11, chr2:10-20, chr2:30-40, chr2:50-60, chr3:10-11</td>
</tr>
</tbody>
</table>
<p>If you walker is implemented in a way that Queue should not divide up your data you should explicitly set the <code>@PartitionBy(PartitionType.NONE)</code>. If your walker can theoretically be run per genome location specify <code>@PartitionBy(PartitionType.LOCUS)</code>.</p>
<pre><code class="pre_md">@PartitionBy(PartitionType.LOCUS)
public class ExampleWalker extends LocusWalker&lt;Integer, Integer&gt; {
...</code class="pre_md"></pre>
<h4>Specifying how to join outputs</h4>
<p>Queue will join the standard walker outputs.</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Output type</th>
<th style="text-align: left;">Default gatherer implementation</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">SAMFileWriter</td>
<td style="text-align: left;">The BAM files are joined together using Picard's MergeSamFiles.</td>
</tr>
<tr>
<td style="text-align: left;">VCFWriter</td>
<td style="text-align: left;">The VCF files are joined together using the GATK CombineVariants.</td>
</tr>
<tr>
<td style="text-align: left;">PrintStream</td>
<td style="text-align: left;">The first two files are scanned for a common header. The header is written once into the output, and then each file is appended to the output, skipping past with the header lines.</td>
</tr>
</tbody>
</table>
<p>If your PrintStream is not a simple text file that can be concatenated together, you must implement a <code>Gatherer</code>. Extend your custom Gatherer from the abstract base class and implement the <code>gather()</code> method.</p>
<pre><code class="pre_md">package org.broadinstitute.sting.commandline;

import java.io.File;
import java.util.List;

/**
 * Combines a list of files into a single output.
 */
public abstract class Gatherer {
    /**
     * Gathers a list of files into a single output.
     * @param inputs Files to combine.
     * @param output Path to output file.
     */
    public abstract void gather(List&lt;File&gt; inputs, File output);

    /**
     * Returns true if the caller should wait for the input files to propagate over NFS before running gather().
     */
    public boolean waitForInputs() { return true; }
}</code class="pre_md"></pre>
<p>Specify your gatherer using the <code>@Gather()</code> annotation by your <code>@Output</code>.</p>
<pre><code class="pre_md">@Output
@Gather(MyGatherer.class)
public PrintStream out;</code class="pre_md"></pre>
<p>Queue will run your custom gatherer to join the intermediate outputs together.</p>
<h3>3. Using GATK walkers in Queue</h3>
<h4>Queue GATK Extensions</h4>
<p>Running 'ant queue' builds a set of Queue extensions for the GATK-Engine. Every GATK walker and command line program in the compiled <code>GenomeAnalysisTK.jar</code> a Queue compatible wrapper is generated.</p>
<p>The extensions can be imported via <code>import org.broadinstitute.sting.queue.extensions.gatk._</code></p>
<pre><code class="pre_md">import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class MyQscript extends QScript {
...</code class="pre_md"></pre>
<p>Note that the generated GATK extensions will automatically handle shell-escaping of all values assigned to the various Walker parameters, so you can rest assured that all of your values will be taken literally by the shell. Do <strong>not</strong> attempt to escape values yourself -- ie.,</p>
<p>Do this: </p>
<pre><code class="pre_md">filterSNPs.filterExpression = List("QD&lt;2.0", "MQ&lt;40.0", "HaplotypeScore&gt;13.0")</code class="pre_md"></pre>
<p>NOT this: </p>
<pre><code class="pre_md">filterSNPs.filterExpression = List("\"QD&lt;2.0\"", "\"MQ&lt;40.0\"", "\"HaplotypeScore&gt;13.0\"")</code class="pre_md"></pre>
<h4>Listing variables</h4>
<p>In addition to the GATK documentation on this wiki you can also find the full list of arguments for each walker extension in a variety of ways.</p>
<p>The source code for the extensions is generated during <code>ant queue</code> and placed in this directory:</p>
<pre><code class="pre_md">build/queue-extensions/src</code class="pre_md"></pre>
<p>When properly configured an IDE can provide command completion of the walker extensions. See <a href="http://gatkforums.broadinstitute.org/discussion/1285/parallelism-with-the-gatk#latest">Queue with IntelliJ IDEA</a> for our recommended settings.</p>
<p>If you do not have access to an IDE you can still find the names of the generated variables using the command line. The generated variable names on each extension are based off of the <code>fullName</code> of the Walker argument. To see the built in documentation for each Walker, run the GATK with:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T &lt;walker name&gt; -help</code class="pre_md"></pre>
<p>Once the import statement is specified you can add() instances of gatk extensions in your QScript's script() method.</p>
<h4>Setting variables</h4>
<p>If the GATK walker input allows more than one of a value you should specify the values as a <code>List()</code>.</p>
<pre><code class="pre_md">  def script() {
    val snps = new UnifiedGenotyper
    snps.reference_file = new File("testdata/exampleFASTA.fasta")
    snps.input_file = List(new File("testdata/exampleBAM.bam"))
    snps.out = new File("snps.vcf")
    add(snps)
  }</code class="pre_md"></pre>
<p>Although it may be harder for others trying to read your QScript, for each of the long name arguments the extensions contain aliases to their short names as well.</p>
<pre><code class="pre_md">  def script() {
    val snps = new UnifiedGenotyper
    snps.R = new File("testdata/exampleFASTA.fasta")
    snps.I = List(new File("testdata/exampleBAM.bam"))
    snps.out = new File("snps.vcf")
    add(snps)
  }</code class="pre_md"></pre>
<p>Here are a few more examples using various list assignment operators.</p>
<pre><code class="pre_md">  def script() {
    val countCovariates = new CountCovariates

    // Append to list using item appender :+
    countCovariates.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)

    // Append to list using collection appender ++
    countCovariates.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")

    // Assign list using plain old object assignment
    countCovariates.input_file = List(inBam)

    // The following is not a list, so just assigning one file to another
    countCovariates.recal_file = outRecalFile

    add(countCovariates)
  }</code class="pre_md"></pre>
<h4>Specifying an alternate GATK jar</h4>
<p>By default Queue runs the GATK from the current classpath. This works best since the extensions are generated and compiled at time same time the GATK is compiled via <code>ant queue</code>.</p>
<p>If you need to swap in a different version of the GATK you may not be able to use the generated extensions. <strong>The alternate GATK jar must have the same command line arguments as the GATK  compiled with Queue.</strong> Otherwise the arguments will not match and you will get an error when Queue attempts to run the alternate GATK jar. In this case you will have to create your own custom <code>CommandLineFunction</code> for your analysis.</p>
<pre><code class="pre_md">  def script {
    val snps = new UnifiedGenotyper
    snps.jarFile = new File("myPatchedGATK.jar")
    snps.reference_file = new File("testdata/exampleFASTA.fasta")
    snps.input_file = List(new File("testdata/exampleBAM.bam"))
    snps.out = new File("snps.vcf")
    add(snps)
  }</code class="pre_md"></pre>
<h4>GATK scatter/gather</h4>
<p>Queue currently allows QScript authors to explicitly invoke scatter/gather on GATK walkers by setting the scatter count on a function.</p>
<pre><code class="pre_md">  def script {
    val snps = new UnifiedGenotyper
    snps.reference_file = new File("testdata/exampleFASTA.fasta")
    snps.input_file = List(new File("testdata/exampleBAM.bam"))
    snps.out = new File("snps.vcf")
    snps.scatterCount = 20
    add(snps)
  }</code class="pre_md"></pre>
<p>This will run the UnifiedGenotyper up to 20 ways parallel and then will merge the partial VCFs back into the single <code>snps.vcf</code>.</p>
<h4>Additional caveat</h4>
<p>Some walkers are still being updated to support Queue fully.  For example they may not have defined the <code>@Input</code> and <code>@Output</code> and thus Queue is unable to correctly track their dependencies, or a custom <code>Gatherer</code> may not be implemented yet.</p>
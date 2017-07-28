## Overview of Queue

http://gatkforums.broadinstitute.org/gatk/discussion/1306/overview-of-queue

<h3>1. Introduction</h3>
<p>GATK-Queue is command-line scripting framework for defining multi-stage genomic analysis pipelines combined with an execution manager that runs those pipelines from end-to-end. Often processing genome data includes several steps to produces outputs, for example our BAM to VCF calling pipeline include among other things:</p>
<ul>
<li>Local realignment around indels</li>
<li>Emitting raw SNP calls</li>
<li>Emitting indels</li>
<li>Masking the SNPs at indels</li>
<li>Annotating SNPs using chip data</li>
<li>Labeling suspicious calls based on filters</li>
<li>Creating a summary report with statistics</li>
</ul>
<p>Running these tools one by one in series may often take weeks for processing, or would require custom scripting to try and optimize using parallel resources.</p>
<p>With a Queue script users can semantically define the multiple steps of the pipeline and then hand off the logistics of running the pipeline to completion. Queue runs independent jobs in parallel, handles transient errors, and uses various techniques such as running multiple copies of the same program on different portions of the genome to produce outputs faster.</p>
<hr />
<h3>2. Obtaining Queue</h3>
<p>You have two options: download the binary distribution (prepackaged, ready to run program) or build it from source.</p>
<h4>- Download the binary</h4>
<p>This is obviously the easiest way to go. Links are on the <a href="http://www.broadinstitute.org/gatk/download">Downloads</a> page. Just get the Queue package; no need to get the GATK package separately as GATK is bundled in with Queue.</p>
<h4>- Building Queue from source</h4>
<p>Briefly, here's what you need to know/do:</p>
<p>Queue is part of the GATK repository. Download the source from the <a href="https://github.com/broadgsa/gatk/">public repository</a> on Github. Run the following command:</p>
<pre><code class="pre_md">git clone https://github.com/broadgsa/gatk.git</code class="pre_md"></pre>
<p><strong>IMPORTANT NOTE:</strong> These instructions refer to the MIT-licensed version of the GATK+Queue source code. With that version, you will be able to build Queue itself, as well as the public portion of the GATK (the core framework), but that will not include the GATK analysis tools. If you want to use Queue to pipeline the GATK analysis tools, you need to clone the <a href="https://github.com/broadgsa/gatk-protected/">'protected' repository</a>. Please note however that part of the source code in that repository (the 'protected' module) is under a different license which excludes for-profit use, modification and redistribution. </p>
<p>Move to the git root directory and use maven to build the source.</p>
<pre><code class="pre_md">mvn clean verify</code class="pre_md"></pre>
<p>All dependencies will be managed by Maven as needed.</p>
<p>See <a href="http://www.broadinstitute.org/gatk/guide/article?id=1287">this article</a> on how to test your installation of Queue.</p>
<hr />
<h3>3. Running Queue</h3>
<p>See <a href="http://www.broadinstitute.org/gatk/guide/article?id=1288">this article</a> on running Queue for the first time for full details.</p>
<p>Queue arguments can be listed by running with <code>--help</code></p>
<pre><code class="pre_md">java -jar dist/Queue.jar --help</code class="pre_md"></pre>
<p>To list the arguments required by a <a href="http://www.broadinstitute.org/gatk/guide/article?id=1307">QScript</a>, add the script with <code>-S</code> and run with <code>--help</code>.</p>
<pre><code class="pre_md">java -jar dist/Queue.jar -S script.scala --help</code class="pre_md"></pre>
<p>Note that by default queue runs in a &quot;dry&quot; mode, as explained in the link above. After verifying the generated commands execute the pipeline by adding <code>-run</code>.</p>
<p>See <a href="http://www.broadinstitute.org/gatk/guide/article?id=1311">QFunction and Command Line Options</a> for more info on adjusting Queue options.</p>
<h3>4. QScripts</h3>
<h4>General Information</h4>
<p>Queue pipelines are written as Scala 2.8 files with a bit of syntactic sugar, called QScripts.</p>
<p>Every QScript includes the following steps:</p>
<ul>
<li>New instances of CommandLineFunctions are created</li>
<li>Input and output arguments are specified on each function</li>
<li>The function is added with <code>add()</code> to Queue for dispatch and monitoring</li>
</ul>
<p>The basic command-line to run the Queue pipelines on the command line is </p>
<pre><code class="pre_md">java -jar Queue.jar -S &lt;script&gt;.scala</code class="pre_md"></pre>
<p>See the main article <a href="http://www.broadinstitute.org/gatk/guide/article?id=1307">Queue QScripts</a> for more info on QScripts.</p>
<h4>Supported QScripts</h4>
<p>Most QScripts are analysis pipelines that are custom-built for specific projects, and we currently do not offer any QScripts as supported analysis tools. However, we do provide some example scripts that you can use as basis to write your own QScripts (see below).</p>
<h4>Example QScripts</h4>
<p>The latest version of the example files are available in the Sting github repository under <a href="https://github.com/broadgsa/gatk/tree/master/public/gatk-queue-extensions-public/src/main/qscripts/org/broadinstitute/gatk/queue/qscripts/examples">public/scala/qscript/examples</a></p>
<hr />
<h3>5. Visualization and Queue</h3>
<h4>QJobReport</h4>
<p>Queue automatically generates <a href="http://www.broadinstitute.org/gatk/guide/article?id=1244">GATKReport</a>-formatted runtime information about executed jobs. See <a href="https://www.dropbox.com/s/jrgba4qojkplk96/QJobReport.pdf?dl=0">this presentation</a> for a general introduction to QJobReport.</p>
<p>Note that Queue attempts to generate a standard visualization using an R script in the GATK <code>public/R</code> repository.  You must provide a path to this location if you want the script to run automatically.  Additionally the script requires the <code>gsalib</code> to be installed on the machine, which is typically done by providing its path in your <code>.Rprofile</code> file:</p>
<pre><code class="pre_md">bm8da-dbe ~/Desktop/broadLocal/GATK/unstable % cat ~/.Rprofile
.libPaths("/Users/depristo/Desktop/broadLocal/GATK/unstable/public/R/")</code class="pre_md"></pre>
<p>Note that gsalib is available from the CRAN repository so you can install it with the canonical R package install command.</p>
<h4>Caveats</h4>
<ul>
<li>
<p>The system only provides information about commands that have just run.  Resuming from a partially completed job will only show the information for the jobs that just ran, and not for any of the completed commands.  This is due to a structural limitation in Queue, and will be fixed when the Queue infrastructure improves</p>
</li>
<li>This feature only works for command line and LSF execution models.  SGE should be easy to add for a motivated individual but we cannot test this capabilities here at the Broad.  Please send us a patch if you do extend Queue to support SGE.</li>
</ul>
<h4>DOT visualization of Pipelines</h4>
<p>Queue emits a <code>queue.dot</code> file to help visualize your commands.  You can open this file in programs like DOT, OmniGraffle, etc to view your pipelines.  By default the system will print out your LSF command lines, but this can be too much in a complex pipeline.  </p>
<p>To clarify your pipeline, override the <code>dotString()</code> function:</p>
<pre><code class="pre_md">class CountCovariates(bamIn: File, recalDataIn: File, args: String = "") extends GatkFunction {
    @Input(doc="foo") var bam = bamIn
    @Input(doc="foo") var bamIndex = bai(bamIn)
    @Output(doc="foo") var recalData = recalDataIn
    memoryLimit = Some(4)
    override def dotString = "CountCovariates: %s [args %s]".format(bamIn.getName, args)
    def commandLine = gatkCommandLine("CountCovariates") + args + " -l INFO -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -I %s --max_reads_at_locus 20000 -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile %s".format(bam, recalData)
}</code class="pre_md"></pre>
<p>Here we only see <code>CountCovariates my.bam [-OQ]</code>, for example, in the <code>dot</code> file.  The base quality score recalibration pipeline, as visualized by DOT, can be viewed here: </p>
<h3>6. Further reading</h3>
<ul>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1288">How to run Queue for the first time</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1309">Queue with IntelliJ IDEA</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1307">Queue QScripts</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1311">QFunction and Command Line Options</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1312">Queue CommandLineFunctions</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1310">Pipelining the GATK using Queue</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1313">Queue with Grid Engine</a></li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=1314">Queue Frequently Asked Questions</a></li>
</ul>
## Writing walkers

http://gatkforums.broadinstitute.org/gatk/discussion/1302/writing-walkers

<h3>1. Introduction</h3>
<p>The core concept behind GATK tools is the walker, a class that implements the three core operations: <strong>filtering</strong>, <strong>mapping</strong>, and <strong>reducing</strong>.</p>
<ul>
<li>
<p><strong>filter</strong>
Reduces the size of the dataset by applying a predicate.</p>
</li>
<li>
<p><strong>map</strong>
Applies a function to each individual element in a dataset, effectively <em>mapping</em> it to a new element.</p>
</li>
<li><strong>reduce</strong>
Inductively combines the elements of a list. The base case is supplied by the <code>reduceInit()</code> function, and the inductive step is performed by the <code>reduce()</code> function.</li>
</ul>
<p>Users of the GATK will provide a walker to run their analyses. The engine will produce a result by first filtering the dataset, running a map operation, and finally reducing the map operation to a single result.</p>
<h3>2. Creating a Walker</h3>
<p>To be usable by the GATK, the walker must satisfy the following properties:</p>
<ul>
<li>
<p>It must subclass one of the basic walkers in the <code>org.broadinstitute.sting.gatk.walkers</code> package, usually ReadWalker or LociWalker.</p>
<ul>
<li>
<p>Locus walkers present all the reads, reference bases, and reference-ordered data that overlap a single base in the reference.  Locus walkers are best used for analyses that look at each locus independently, such as genotyping.</p>
</li>
<li>
<p>Read walkers present only one read at a time, as well as the reference bases and reference-ordered data that overlap that read.</p>
</li>
<li>Besides read walkers and locus walkers, the GATK features several other data access patterns, described <a href="http://www.broadinstitute.org/gatk/guide/article?id=1351">here</a>.</li>
</ul>
</li>
<li>The compiled class or jar must be on the current classpath.  The Java classpath can be controlled using either the <code>$CLASSPATH</code> environment variable or the JVM's <code>-cp</code> option.</li>
</ul>
<h3>3. Examples</h3>
<p>The best way to get started with the GATK is to explore the walkers we've written.  Here are the best walkers to look at when getting started:</p>
<ul>
<li>
<p>CountLoci </p>
<p>It is the simplest locus walker in our codebase. It counts the number of loci walked over in a single run of the GATK.</p>
</li>
</ul>
<p><code>$STING_HOME/java/src/org/broadinstitute/sting/gatk/walkers/qc/CountLociWalker.java</code></p>
<ul>
<li>
<p>CountReads </p>
<p>It is the simplest read walker in our codebase. It counts the number of reads walked over in a single run of the GATK.</p>
</li>
</ul>
<p><code>$STING_HOME/java/src/org/broadinstitute/sting/gatk/walkers/qc/CountReadsWalker.java</code></p>
<ul>
<li>
<p>GATKPaperGenotyper </p>
<p>This is a more sophisticated example, taken from our recent paper in Genome Research (and using our ReadBackedPileup to select and filter reads). It is an extremely basic Bayesian genotyper that demonstrates how to output data to a stream and execute simple base operations.</p>
</li>
</ul>
<p><code>$STING_HOME/java/src/org/broadinstitute/sting/gatk/examples/papergenotyper/GATKPaperGenotyper.java</code> </p>
<p><strong>Please note that the walker above is NOT the UnifiedGenotyper.  While conceptually similar to the UnifiedGenotyper, the GATKPaperGenotyper uses a much simpler calling model for increased clarity and readability.</strong></p>
<h3>4. External walkers and the 'external' directory</h3>
<p>The GATK can absorb external walkers placed in a directory of your choosing.  By default, that directory is called 'external' and is relative to the Sting git root directory (for example, <code>~/src/Sting/external</code>).  However, you can choose to place that directory anywhere on the filesystem and specify its complete path using the ant <code>external.dir</code> property. </p>
<pre><code class="pre_md">ant -Dexternal.dir=~/src/external</code class="pre_md"></pre>
<p>The GATK will check each directory under the external directory (but not the external directory itself!) for small build scripts.  These build scripts must contain at least a <code>compile</code> target that compiles your walker and places the resulting class file into the GATK's class file output directory.  The following is a sample compile target:</p>
<pre><code class="pre_md">&lt;target name="compile" depends="init"&gt;
    &lt;javac srcdir="." destdir="${build.dir}" classpath="${gatk.classpath}" /&gt;
&lt;/target&gt;</code class="pre_md"></pre>
<p>As a convenience, the <code>build.dir</code> ant property will be predefined to be the GATK's class file output directory and the <code>gatk.classpath</code> property will be predefined to be the GATK's core classpath.  Once this structure is defined, any invocation of the ant build scripts will build the contents of the external directory as well as the GATK itself.</p>
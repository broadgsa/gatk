## Writing walkers in Scala

http://gatkforums.broadinstitute.org/gatk/discussion/1354/writing-walkers-in-scala

<h2>1. Install scala somewhere</h2>
<p>At the Broad, we typically put it somewhere like this: </p>
<pre><code class="pre_md">/home/radon01/depristo/work/local/scala-2.7.5.final</code class="pre_md"></pre>
<p>Next, create a symlink from this directory to <code>trunk/scala/installation</code>:</p>
<pre><code class="pre_md">ln -s /home/radon01/depristo/work/local/scala-2.7.5.final trunk/scala/installation</code class="pre_md"></pre>
<h2>2. Setting up your path</h2>
<p>Right now the only way to get scala walkers into the GATK is by explicitly setting your <code>CLASSPATH</code> in your <code>.my.cshrc</code> file:</p>
<pre><code class="pre_md">setenv CLASSPATH /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/FourBaseRecaller.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/Playground.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/StingUtils.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/bcel-5.2.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/colt-1.2.0.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/google-collections-0.9.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/javassist-3.7.ga.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/junit-4.4.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/log4j-1.2.15.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/picard-1.02.63.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/picard-private-875.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/reflections-0.9.2.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/sam-1.01.63.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/simple-xml-2.0.4.jar:/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/GATKScala.jar:/humgen/gsa-scr1/depristo/local/scala-2.7.5.final/lib/scala-library.jar</code class="pre_md"></pre>
<p>Really this needs to be manually updated whenever any of the libraries are updated.  If you see this error:</p>
<pre><code class="pre_md">Caused by: java.lang.RuntimeException: java.util.zip.ZipException: error in opening zip file
        at org.reflections.util.VirtualFile.iterable(VirtualFile.java:79)
        at org.reflections.util.VirtualFile$5.transform(VirtualFile.java:169)
        at org.reflections.util.VirtualFile$5.transform(VirtualFile.java:167)
        at org.reflections.util.FluentIterable$3.transform(FluentIterable.java:43)
        at org.reflections.util.FluentIterable$3.transform(FluentIterable.java:41)
        at org.reflections.util.FluentIterable$ForkIterator.computeNext(FluentIterable.java:81)
        at com.google.common.collect.AbstractIterator.tryToComputeNext(AbstractIterator.java:132)
        at com.google.common.collect.AbstractIterator.hasNext(AbstractIterator.java:127)
        at org.reflections.util.FluentIterable$FilterIterator.computeNext(FluentIterable.java:102)
        at com.google.common.collect.AbstractIterator.tryToComputeNext(AbstractIterator.java:132)
        at com.google.common.collect.AbstractIterator.hasNext(AbstractIterator.java:127)
        at org.reflections.util.FluentIterable$TransformIterator.computeNext(FluentIterable.java:124)
        at com.google.common.collect.AbstractIterator.tryToComputeNext(AbstractIterator.java:132)
        at com.google.common.collect.AbstractIterator.hasNext(AbstractIterator.java:127)
        at org.reflections.Reflections.scan(Reflections.java:69)
        at org.reflections.Reflections.&lt;init&gt;(Reflections.java:47)
        at org.broadinstitute.sting.utils.PackageUtils.&lt;clinit&gt;(PackageUtils.java:23)</code class="pre_md"></pre>
<p>It's because the libraries aren't updated.  Basically just do an <code>ls</code> of your <code>trunk/dist</code> directory after the GATK has been build, make this your classpath as above, and tack on:</p>
<pre><code class="pre_md">/humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/GATKScala.jar:/humgen/gsa-scr1/depristo/local/scala-2.7.5.final/lib/scala-library.jar</code class="pre_md"></pre>
<p>A command that almost works (but you'll need to replace the spaces with colons) is:</p>
<pre><code class="pre_md">#setenv CLASSPATH $CLASSPATH `ls /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/*.jar` /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/GATKScala.jar:/humgen/gsa-scr1/depristo/local/scala-2.7.5.final/lib/scala-library.jar</code class="pre_md"></pre>
<h2>3. Building scala code</h2>
<p>All of the Scala source code lives in <code>scala/src</code>, which you build using <code>ant scala</code></p>
<p>There are already some example Scala walkers in <code>scala/src</code>, so doing a standard checkout, installing scala, settting up your environment, should allow you to run something like:</p>
<pre><code class="pre_md">gsa2 ~/dev/GenomeAnalysisTK/trunk &gt; ant scala
Buildfile: build.xml

init.scala:

scala:
     [echo] Sting: Compiling scala!
   [scalac] Compiling 2 source files to /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/scala/classes
   [scalac] warning: there were deprecation warnings; re-run with -deprecation for details
   [scalac] one warning found
   [scalac] Compile suceeded with 1 warning; see the compiler output for details.
   [delete] Deleting: /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/GATKScala.jar
      [jar] Building jar: /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/dist/GATKScala.jar</code class="pre_md"></pre>
<h2>4. Invoking a scala walker</h2>
<p>Until we can include Scala walkers along with the main GATK jar (avoiding the classpath issue too) you have to invoke your scala walkers using this syntax:</p>
<pre><code class="pre_md">java -Xmx2048m org.broadinstitute.sting.gatk.CommandLineGATK -T BaseTransitionTableCalculator -R /broad/1KG/reference/human_b36_both.fasta -I /broad/1KG/DCC_merged/freeze5/NA12878.pilot2.SLX.bam -l INFO -L 1:1-100</code class="pre_md"></pre>
<p>Here, the <code>BaseTransitionTableCalculator</code> walker is written in Scala and being loaded into the system by the GATK walker manager.  Otherwise everything looks like a normal GATK module.</p>
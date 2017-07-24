## How to access the picard and htsjdk repository (containing samtools-jdk, tribble, and variant)

http://gatkforums.broadinstitute.org/gatk/discussion/2194/how-to-access-the-picard-and-htsjdk-repository-containing-samtools-jdk-tribble-and-variant

<p>The picard repository on github contains all picard public tools. Libraries live under the htsjdk, which includes the samtools-jdk, tribble, and variant packages (which includes VariantContext and associated classes as well as the VCF/BCF codecs).</p>
<p>If you just need to check out the sources and don't need to make any commits into the picard repository, the command is:</p>
<pre><code class="pre_md">git clone https://github.com/broadinstitute/picard.git</code class="pre_md"></pre>
<p>Then within the picard directory, clone the htsjdk.</p>
<pre><code class="pre_md">cd picard
git clone https://github.com/samtools/htsjdk.git</code class="pre_md"></pre>
<p>Then you can attach the <code>picard/src/java</code> and <code>picard/htsjdk/src/java</code> directories in IntelliJ as a source directory (File -&gt; Project Structure -&gt; Libraries -&gt; Click the plus sign -&gt; &quot;Attach Files or Directories&quot; in the latest IntelliJ).</p>
<p>To build picard and the htsjdk all at once, type <code>ant</code> from within the picard directory. To run tests, type <code>ant test</code></p>
<p>If you do need to make commits into the picard repository, first you'll need to create a github account, fork picard or htsjdk, make your changes, and then issue a pull request. For more info on pull requests, see: <a href="https://help.github.com/articles/using-pull-requests">https://help.github.com/articles/using-pull-requests</a></p>
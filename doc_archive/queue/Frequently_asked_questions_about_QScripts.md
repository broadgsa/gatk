## Frequently asked questions about QScripts

http://gatkforums.broadinstitute.org/gatk/discussion/1314/frequently-asked-questions-about-qscripts

<h3>1. Many of my GATK functions are setup with the same Reference, Intervals, etc. Is there a quick way to reuse these values for the different analyses in my pipeline?</h3>
<p>Yes.</p>
<ul>
<li>Create a trait that extends from CommandLineGATK.</li>
<li>In the trait, copy common values from your qscript.</li>
<li>Mix the trait into instances of your classes.</li>
</ul>
<p>For more information, see the <code>ExampleUnifiedGenotyper.scala</code> or examples of using Scala's traits/mixins illustrated in the <a href="http://gatkforums.broadinstitute.org/discussion/1307/queue-pipeline-scripts-qscripts">QScripts documentation</a>.</p>
<h3>2. How do I accept a list of arguments to my QScript?</h3>
<p>In your QScript, define a <code>var</code> list and annotate it with <code>@Argument</code>. Initialize the value to <code>Nil</code>.</p>
<pre><code class="pre_md">@Argument(doc="filter names", shortName="filter")
var filterNames: List[String] = Nil</code class="pre_md"></pre>
<p>On the command line specify the arguments by repeating the argument name.</p>
<pre><code class="pre_md">-filter filter1 -filter filter2 -filter filter3</code class="pre_md"></pre>
<p>Then once your QScript is run, the command line arguments will be available for use in the QScript's <code>script</code> method.</p>
<pre><code class="pre_md">  def script {
     var myCommand = new MyFunction
     myCommand.filters = this.filterNames
  }</code class="pre_md"></pre>
<p>For a full example of command line arguments see the <a href="http://gatkforums.broadinstitute.org/discussion/1307/queue-pipeline-scripts-qscripts">QScripts documentation</a>.</p>
<h3>3. What is the best way to run a utility method at the right time?</h3>
<p>Wrap the utility with an <code>InProcessFunction</code>. If your functionality is reusable code you should add it to <code>Sting Utils</code> with <code>Unit Tests</code> and then invoke your new function from your <code>InProcessFunction</code>. Computationally or memory intensive functions should NOT be implemented as <code>InProcessFunctions</code>, and should be wrapped in <a href="http://gatkforums.broadinstitute.org/discussion/1312/queue-commandlinefunctions">Queue CommandLineFunctions</a> instead.</p>
<pre><code class="pre_md">    class MySplitter extends InProcessFunction {
      @Input(doc="inputs")
      var in: File = _

      @Output(doc="outputs")
      var out: List[File] = Nil

      def run {
         StingUtilityMethod.quickSplitFile(in, out)
      }
    }

    var splitter = new MySplitter
    splitter.in = new File("input.txt")
    splitter.out = List(new File("out1.txt"), new File("out2.txt"))
    add(splitter)</code class="pre_md"></pre>
<p>See <a href="http://gatkforums.broadinstitute.org/discussion/1312/queue-commandlinefunctions">Queue CommandLineFunctions</a> for more information on how <code>@Input</code> and <code>@Output</code> are used.</p>
<h3>4. What is the best way to write a list of files?</h3>
<p>Create an instance of a <code>ListWriterFunction</code> and add it in your script method. </p>
<pre><code class="pre_md">import org.broadinstitute.sting.queue.function.ListWriterFunction

val writeBamList = new ListWriterFunction
writeBamList.inputFiles = bamFiles
writeBamList.listFile = new File("myBams.list")
add(writeBamList)</code class="pre_md"></pre>
<h3>5. How do I add optional debug output to my QScript?</h3>
<p>Queue contains a trait mixin you can use to add Log4J support to your classes.</p>
<p>Add the import for the trait <code>Logging</code> to your QScript.</p>
<pre><code class="pre_md">import org.broadinstitute.sting.queue.util.Logging</code class="pre_md"></pre>
<p>Mixin the trait to your class.</p>
<pre><code class="pre_md">class MyScript extends Logging {
...</code class="pre_md"></pre>
<p>Then use the mixed in <code>logger</code> to write debug output when the user specifies <code>-l DEBUG</code>.</p>
<pre><code class="pre_md">logger.debug("This will only be displayed when debugging is enabled.")</code class="pre_md"></pre>
<h3>6. I updated Queue and now I'm getting java.lang.NoClassDefFoundError / java.lang.AbstractMethodError</h3>
<p>Try <code>ant clean</code>.</p>
<p>Queue relies on a lot of Scala traits / mixins. These dependencies are not always picked up by the scala/java compilers leading to partially implemented classes. If that doesn't work please let us know in the <a href="http://gatkforums.broadinstitute.org/">forum</a>.</p>
<h3>7. Do I need to create directories in my QScript?</h3>
<p>No. QScript will create all parent directories for outputs.</p>
<h3>8. How do I specify the -W 240 for the LSF hour queue at the Broad?</h3>
<p>Queue's LSF dispatcher automatically looks up and sets the maximum runtime for whichever LSF queue is specified. If you set your <code>-jobQueue/.jobQueue</code> to <code>hour</code> then you should see something like this under <code>bjobs -l</code>:</p>
<pre><code class="pre_md">RUNLIMIT
240.0 min of gsa3</code class="pre_md"></pre>
<h3>9. Can I run Queue with GridEngine?</h3>
<p>Queue GridEngine functionality is community supported. See here for full details: <a href="http://gatkforums.broadinstitute.org/discussion/1313/queue-with-grid-engine">Queue with Grid Engine</a>.</p>
<h3>10. How do I pass advanced java arguments to my GATK commands, such as remote debugging?</h3>
<p>The easiest way to do this at the moment is to mixin a trait.</p>
<p>First define a trait which adds your java options:</p>
<pre><code class="pre_md">  trait RemoteDebugging extends JavaCommandLineFunction {
    override def javaOpts = super.javaOpts + " -Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=5005"
  }</code class="pre_md"></pre>
<p>Then mix in the trait to your walker and otherwise run it as normal:</p>
<pre><code class="pre_md">  val printReadsDebug = new PrintReads with RemoteDebugging
  printReadsDebug.reference_sequence = "my.fasta"
  // continue setting up your walker...
  add(printReadsDebug)</code class="pre_md"></pre>
<h3>11. Why does Queue log &quot;Running jobs. ... Done.&quot; but doesn't actually run anything?</h3>
<p>If you see something like the following, it means that Queue believes that it previously successfully generated all of the outputs.</p>
<pre><code class="pre_md">INFO 16:25:55,049 QCommandLine - Scripting ExampleUnifiedGenotyper 
INFO 16:25:55,140 QCommandLine - Added 4 functions 
INFO 16:25:55,140 QGraph - Generating graph. 
INFO 16:25:55,164 QGraph - Generating scatter gather jobs. 
INFO 16:25:55,714 QGraph - Removing original jobs. 
INFO 16:25:55,716 QGraph - Adding scatter gather jobs. 
INFO 16:25:55,779 QGraph - Regenerating graph. 
INFO 16:25:55,790 QGraph - Running jobs. 
INFO 16:25:55,853 QGraph - 0 Pend, 0 Run, 0 Fail, 10 Done 
INFO 16:25:55,902 QCommandLine - Done </code class="pre_md"></pre>
<p>Queue will not re-run the job if a <code>.done</code> file is found for the all the outputs, <em>e.g.</em>: <code>/path/to/.output.file.done</code>. You can either remove the specific <code>.done</code> files yourself, or use the <code>-startFromScratch</code> command line option.</p>
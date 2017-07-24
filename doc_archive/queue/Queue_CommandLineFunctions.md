## Queue CommandLineFunctions

http://gatkforums.broadinstitute.org/gatk/discussion/1312/queue-commandlinefunctions

<h3>1. Basic QScript run rules</h3>
<ul>
<li>In the <code>script</code> method, a QScript will add one or more <code>CommandLineFunction</code>s.</li>
<li>Queue tracks dependencies between functions via variables annotated with <code>@Input</code> and <code>@Output</code>.</li>
<li>Queue will run functions based on the dependencies between them, so if the <code>@Input</code> of <code>CommandLineFunction</code> <code>A</code> depends on the <code>@Output</code> of <code>ComandLineFunction</code> <code>B</code>, <code>A</code> will wait for <code>B</code> to finish before it starts running.</li>
</ul>
<h3>2. Command Line</h3>
<p>Each CommandLineFunction must define the actual command line to run as follows.</p>
<pre><code class="pre_md">class MyCommandLine extends CommandLineFunction {
  def commandLine = "myScript.sh hello world"
}</code class="pre_md"></pre>
<h4>Constructing a Command Line Manually</h4>
<p>If you're writing a one-off CommandLineFunction that is not destined for use
by other QScripts, it's often easiest to construct the command line directly
rather than through the API methods provided in the CommandLineFunction class. </p>
<p>For example:</p>
<pre><code class="pre_md">def commandLine = "cat %s | grep -v \"#\" &gt; %s".format(files, out)</code class="pre_md"></pre>
<h4>Constructing a Command Line using API Methods</h4>
<p>If you're writing a CommandLineFunction that will become part of Queue and/or
will be used by other QScripts, however, our best practice recommendation is
to construct your command line <em>only</em> using the methods provided in the
CommandLineFunction class: <code>required()</code>, <code>optional()</code>, <code>conditional()</code>, and <code>repeat()</code></p>
<p>The reason for this is that these methods automatically escape the values you
give them so that they'll be interpreted literally within the shell scripts
Queue generates to run your command, and they also manage whitespace separation of command-line tokens for you. This prevents (for example) a value like <code>MQ &gt; 10</code> from being interpreted as an output redirection by the shell, and avoids issues with values containing embedded spaces. The methods also give you the ability to turn escaping and/or whitespace separation off as needed. An example:</p>
<pre><code class="pre_md">override def commandLine = super.commandLine +
                           required("eff") +
                           conditional(verbose, "-v") +
                           optional("-c", config) +
                           required("-i", "vcf") +
                           required("-o", "vcf") +
                           required(genomeVersion) +
                           required(inVcf) +
                           required("&gt;", escape=false) +  // This will be shell-interpreted as an output redirection
                           required(outVcf)</code class="pre_md"></pre>
<p>The CommandLineFunctions built into Queue, including the CommandLineFunctions
automatically generated for GATK Walkers, are all written using this pattern.
This means that when you configure a GATK Walker or one of the other built-in
CommandLineFunctions in a QScript, you can rely on all of your values being
safely escaped and taken literally when the commands are run, including values
containing characters that would normally be interpreted by the shell such as
<code>MQ &gt; 10</code>.</p>
<p>Below is a brief overview of the API methods available to you in the <code>CommandLineFunction</code> class for safely constructing command lines:</p>
<ul>
<li><code>required()</code> </li>
</ul>
<p>Used for command-line arguments that are always present, <em>e.g.</em>:</p>
<pre><code class="pre_md">required("-f", "filename")                              returns: " '-f' 'filename' "
required("-f", "filename", escape=false)                returns: " -f filename "
required("java")                                        returns: " 'java' "
required("INPUT=", "myBam.bam", spaceSeparated=false)   returns: " 'INPUT=myBam.bam' "</code class="pre_md"></pre>
<ul>
<li><code>optional()</code> </li>
</ul>
<p>Used for command-line arguments that may or may not be present, <em>e.g.</em>:</p>
<pre><code class="pre_md">optional("-f", myVar) behaves like required() if myVar has a value, but returns ""
if myVar is null/Nil/None</code class="pre_md"></pre>
<ul>
<li><code>conditional()</code> </li>
</ul>
<p>Used for command-line arguments that should only be included if some condition is true, <em>e.g.</em>:</p>
<pre><code class="pre_md">conditional(verbose, "-v") returns " '-v' " if verbose is true, otherwise returns ""</code class="pre_md"></pre>
<ul>
<li><code>repeat()</code> </li>
</ul>
<p>Used for command-line arguments that are repeated multiple times on the command line, <em>e.g.</em>:</p>
<pre><code class="pre_md">repeat("-f", List("file1", "file2", "file3")) returns: " '-f' 'file1' '-f' 'file2' '-f' 'file3' "</code class="pre_md"></pre>
<h3>3. Arguments</h3>
<ul>
<li>
<p><code>CommandLineFunction</code> arguments use a similar syntax to arguments.</p>
</li>
<li><code>CommandLineFunction</code> variables are annotated with <code>@Input</code>, <code>@Output</code>, or <code>@Argument</code> annotations.</li>
</ul>
<h4>Input and Output Files</h4>
<p>So that Queue can track the input and output files of a command, <code>CommandLineFunction</code> <code>@Input</code> and <code>@Output</code> must be <code>java.io.File</code> objects.</p>
<pre><code class="pre_md">class MyCommandLine extends CommandLineFunction {
  @Input(doc="input file")
  var inputFile: File = _
  def commandLine = "myScript.sh -fileParam " + inputFile
}</code class="pre_md"></pre>
<h4>FileProvider</h4>
<p><code>CommandLineFunction</code> variables can also provide indirect access to <code>java.io.File</code> inputs and outputs via the <code>FileProvider</code> trait.</p>
<pre><code class="pre_md">class MyCommandLine extends CommandLineFunction {
  @Input(doc="named input file")
  var inputFile: ExampleFileProvider = _
  def commandLine = "myScript.sh " + inputFile
}

// An example FileProvider that stores a 'name' with a 'file'.
class ExampleFileProvider(var name: String, var file: File) extends org.broadinstitute.sting.queue.function.FileProvider {
  override def toString = " -fileName " + name + " -fileParam " + file
}</code class="pre_md"></pre>
<h4>Optional Arguments</h4>
<p>Optional files can be specified via <code>required=false</code>, and can use the <code>CommandLineFunction.optional()</code> utility method, as described above:</p>
<pre><code class="pre_md">class MyCommandLine extends CommandLineFunction {
  @Input(doc="input file", required=false)
  var inputFile: File = _
  // -fileParam will only be added if the QScript sets inputFile on this instance of MyCommandLine
  def commandLine = required("myScript.sh") + optional("-fileParam", inputFile)
}</code class="pre_md"></pre>
<h4>Collections as Arguments</h4>
<p>A <code>List</code> or <code>Set</code> of files can use the <code>CommandLineFunction.repeat()</code> utility method, as described above:</p>
<pre><code class="pre_md">class MyCommandLine extends CommandLineFunction {
  @Input(doc="input file")
  var inputFile: List[File] = Nil // NOTE: Do not set List or Set variables to null!
  // -fileParam will added as many times as the QScript adds the inputFile on this instance of MyCommandLine
  def commandLine = required("myScript.sh") + repeat("-fileParam", inputFile)
}</code class="pre_md"></pre>
<h4>Non-File Arguments</h4>
<p>A command line function can define other required arguments via @Argument.</p>
<pre><code class="pre_md">class MyCommandLine extends CommandLineFunction {
  @Argument(doc="message to display")
  var veryImportantMessage: String = _
  // If the QScript does not specify the required veryImportantMessage, the pipeline will not run.
  def commandLine = required("myScript.sh") + required(veryImportantMessage)
}</code class="pre_md"></pre>
<h3>4. Example: &quot;samtools index&quot;</h3>
<pre><code class="pre_md">class SamToolsIndex extends CommandLineFunction {
  @Input(doc="bam to index") var bamFile: File = _
  @Output(doc="bam index") var baiFile: File = _
  def commandLine = "samtools index %s %s".format(bamFile, baiFile)
)</code class="pre_md"></pre>
<p>Or, using the CommandLineFunction API methods to construct the command line with automatic shell escaping:</p>
<pre><code class="pre_md">class SamToolsIndex extends CommandLineFunction {
  @Input(doc="bam to index") var bamFile: File = _
  @Output(doc="bam index") var baiFile: File = _
  def commandLine = required("samtools") + required("index") + required(bamFile) + required(baiFile)
)</code class="pre_md"></pre>
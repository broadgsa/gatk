## Queue pipeline scripts (QScripts)

http://gatkforums.broadinstitute.org/gatk/discussion/1307/queue-pipeline-scripts-qscripts

<h3>1. Introduction</h3>
<p>Queue pipelines are Scala 2.8 files with a bit of syntactic sugar, called QScripts. Check out the following as references.</p>
<ul>
<li><a href="http://programming-scala.labs.oreilly.com">http://programming-scala.labs.oreilly.com</a></li>
<li><a href="http://www.scala-lang.org/docu/files/ScalaByExample.pdf">http://www.scala-lang.org/docu/files/ScalaByExample.pdf</a></li>
<li><a href="http://davetron5000.github.com/scala-style/index.html">http://davetron5000.github.com/scala-style/index.html</a></li>
</ul>
<p>QScripts are easiest to develop using an Integrated Development Environment. See <a href="http://gatkforums.broadinstitute.org/discussion/1309/queue-with-intellij-idea">Queue with IntelliJ IDEA</a> for our recommended settings.</p>
<p>The following is a basic outline of a QScript:</p>
<pre><code class="pre_md">import org.broadinstitute.sting.queue.QScript
// List other imports here

// Define the overall QScript here.
class MyScript extends QScript {
  // List script arguments here.
  @Input(doc="My QScript inputs")
  var scriptInput: File = _

  // Create and add the functions in the script here.
  def script = {
     var myCL = new MyCommandLine
     myCL.myInput = scriptInput // Example variable input
     myCL.myOutput = new File("/path/to/output") // Example hardcoded output
     add(myCL)
  }

}</code class="pre_md"></pre>
<h3>2. Imports</h3>
<p>Imports can be any scala or java imports in scala syntax.</p>
<pre><code class="pre_md">import java.io.File
import scala.util.Random
import org.favorite.my._
// etc.</code class="pre_md"></pre>
<h3>3. Classes</h3>
<ul>
<li>
<p>To add a <code>CommandLineFunction</code> to a pipeline, a class must be defined that extends <code>QScript</code>.</p>
</li>
<li>
<p>The <code>QScript</code> must define a method <code>script</code>.</p>
</li>
<li>The <code>QScript</code> can define helper methods or variables.</li>
</ul>
<h3>4. Script method</h3>
<p>The body of <code>script</code> should create and add <a href="http://gatkforums.broadinstitute.org/discussion/1312/queue-commandlinefunctions">Queue CommandlineFunctions</a>.</p>
<pre><code class="pre_md">class MyScript extends org.broadinstitute.sting.queue.QScript {
  def script = add(new CommandLineFunction { def commandLine = "echo hello world" })
}</code class="pre_md"></pre>
<h3>5. Command Line Arguments</h3>
<ul>
<li>
<p>A <code>QScript</code> canbe set to read command line arguments by defining variables with <code>@Input</code>, <code>@Output</code>, or <code>@Argument</code> annotations.</p>
</li>
<li>
<p>A command line argument can be a primitive scalar, enum, <code>File</code>, or scala immutable <code>Array</code>, <code>List</code>, <code>Set</code>, or <code>Option</code> of a primitive, enum, or <code>File</code>.</p>
</li>
<li>
<p><code>QScript</code> command line arguments can be marked as optional by setting <code>required=false</code>.</p>
<p>class MyScript extends org.broadinstitute.sting.queue.QScript {
@Input(doc=&quot;example message to echo&quot;)
var message: String = _
def script = add(new CommandLineFunction { def commandLine = &quot;echo &quot; + message })
}</p>
</li>
</ul>
<h3>6. Using and writing CommandLineFunctions</h3>
<h4>Adding existing GATK walkers</h4>
<p>See <a href="http://gatkforums.broadinstitute.org/discussion/1310/pipelining-the-gatk-with-queue">Pipelining the GATK using Queue</a> for more information on the automatically generated Queue wrappers for GATK walkers.</p>
<p>After functions are defined they should be added to the <code>QScript</code> pipeline using <code>add()</code>.</p>
<pre><code class="pre_md">for (vcf &lt;- vcfs) {
  val ve = new VariantEval
  ve.vcfFile = vcf
  ve.evalFile = swapExt(vcf, "vcf", "eval")
  add(ve)
}</code class="pre_md"></pre>
<h4>Defining new CommandLineFunctions</h4>
<ul>
<li>
<p>Queue tracks dependencies between functions via variables annotated with <code>@Input</code> and <code>@Output</code>.</p>
</li>
<li>
<p>Queue will run functions based on the dependencies between them, not based on the order in which they are added in the script! So if the <code>@Input</code> of <code>CommandLineFunction</code> <code>A</code> depends on the <code>@Output</code> of <code>ComandLineFunction</code> <code>B</code>, <code>A</code> will wait for <code>B</code> to finish before it starts running.</p>
</li>
<li>See the main article <a href="http://gatkforums.broadinstitute.org/discussion/1312/queue-commandlinefunctions">Queue CommandLineFunctions</a> for more information.</li>
</ul>
<h3>7. Examples</h3>
<ul>
<li>
<p>The latest version of the example files are available in the Sting git repository under [public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/](<a href="https://github.com/broadgsa/gatk/blob/master/public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/">https://github.com/broadgsa/gatk/blob/master/public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/</a> ).</p>
</li>
<li>To print the list of arguments required by an existing QScript run with <code>-help</code>.</li>
<li>To check if your script has all of the <code>CommandLineFunction</code> variables set correctly, run <em>without</em> <code>-run</code>.</li>
<li>When you are ready to execute the full pipeline, add <code>-run</code>.</li>
</ul>
<h4>Hello World QScript</h4>
<p>The following is a &quot;hello world&quot; example that runs a single command line to <code>echo hello world</code>.</p>
<pre><code class="pre_md">import org.broadinstitute.sting.queue.QScript

class HelloWorld extends QScript {
  def script = {
    add(new CommandLineFunction {
      def commandLine = "echo hello world"
    })
  }
}</code class="pre_md"></pre>
<p>The above file is checked into the Sting git repository under <a href="https://github.com/broadgsa/gatk/blob/master/public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala">HelloWorld.scala</a>. After building Queue from source, the QScript can be run with the following command:</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar dist/Queue.jar -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala -run</code class="pre_md"></pre>
<p>It should produce output similar to:</p>
<pre><code class="pre_md">INFO  16:23:27,825 QScriptManager - Compiling 1 QScript 
INFO  16:23:31,289 QScriptManager - Compilation complete 
INFO  16:23:34,631 HelpFormatter - --------------------------------------------------------- 
INFO  16:23:34,631 HelpFormatter - Program Name: org.broadinstitute.sting.queue.QCommandLine 
INFO  16:23:34,632 HelpFormatter - Program Args: -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala -run  
INFO  16:23:34,632 HelpFormatter - Date/Time: 2011/01/14 16:23:34 
INFO  16:23:34,632 HelpFormatter - --------------------------------------------------------- 
INFO  16:23:34,632 HelpFormatter - --------------------------------------------------------- 
INFO  16:23:34,634 QCommandLine - Scripting HelloWorld 
INFO  16:23:34,651 QCommandLine - Added 1 functions 
INFO  16:23:34,651 QGraph - Generating graph. 
INFO  16:23:34,660 QGraph - Running jobs. 
INFO  16:23:34,689 ShellJobRunner - Starting: echo hello world 
INFO  16:23:34,689 ShellJobRunner - Output written to /Users/kshakir/src/Sting/Q-43031@bmef8-d8e-1.out 
INFO  16:23:34,771 ShellJobRunner - Done: echo hello world 
INFO  16:23:34,773 QGraph - Deleting intermediate files. 
INFO  16:23:34,773 QCommandLine - Done </code class="pre_md"></pre>
<h4>ExampleUnifiedGenotyper.scala</h4>
<p>This example uses automatically generated Queue compatible wrappers for the GATK. See <a href="http://gatkforums.broadinstitute.org/discussion/1310/pipelining-the-gatk-with-queue">Pipelining the GATK using Queue</a> for more info on authoring Queue support into walkers and using walkers in Queue.</p>
<p>The <a href="https://github.com/broadgsa/gatk/blob/master/public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala">ExampleUnifiedGenotyper.scala</a> for running the UnifiedGenotyper followed by VariantFiltration can be found in the examples folder.</p>
<p>To list the command line parameters, including the required parameters, run with <code>-help</code>.</p>
<pre><code class="pre_md">java -jar dist/Queue.jar -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala -help</code class="pre_md"></pre>
<p>The help output should appear similar to this:</p>
<pre><code class="pre_md">INFO  10:26:08,491 QScriptManager - Compiling 1 QScript
INFO  10:26:11,926 QScriptManager - Compilation complete
---------------------------------------------------------
Program Name: org.broadinstitute.sting.queue.QCommandLine
---------------------------------------------------------
---------------------------------------------------------
usage: java -jar Queue.jar -S &lt;script&gt; [-run] [-jobRunner &lt;job_runner&gt;] [-bsub] [-status] [-retry &lt;retry_failed&gt;]
       [-startFromScratch] [-keepIntermediates] [-statusTo &lt;status_email_to&gt;] [-statusFrom &lt;status_email_from&gt;] [-dot
       &lt;dot_graph&gt;] [-expandedDot &lt;expanded_dot_graph&gt;] [-jobPrefix &lt;job_name_prefix&gt;] [-jobProject &lt;job_project&gt;] [-jobQueue
       &lt;job_queue&gt;] [-jobPriority &lt;job_priority&gt;] [-memLimit &lt;default_memory_limit&gt;] [-runDir &lt;run_directory&gt;] [-tempDir
       &lt;temp_directory&gt;] [-jobSGDir &lt;job_scatter_gather_directory&gt;] [-emailHost &lt;emailSmtpHost&gt;] [-emailPort &lt;emailSmtpPort&gt;]
       [-emailTLS] [-emailSSL] [-emailUser &lt;emailUsername&gt;] [-emailPassFile &lt;emailPasswordFile&gt;] [-emailPass &lt;emailPassword&gt;]
       [-l &lt;logging_level&gt;] [-log &lt;log_to_file&gt;] [-quiet] [-debug] [-h] -R &lt;referencefile&gt; -I &lt;bamfile&gt; [-L &lt;intervals&gt;]
       [-filter &lt;filternames&gt;] [-filterExpression &lt;filterexpressions&gt;]

 -S,--script &lt;script&gt;                                                      QScript scala file
 -run,--run_scripts                                                        Run QScripts.  Without this flag set only
                                                                           performs a dry run.
 -jobRunner,--job_runner &lt;job_runner&gt;                                      Use the specified job runner to dispatch
                                                                           command line jobs
 -bsub,--bsub                                                              Equivalent to -jobRunner Lsf706
 -status,--status                                                          Get status of jobs for the qscript
 -retry,--retry_failed &lt;retry_failed&gt;                                      Retry the specified number of times after a
                                                                           command fails.  Defaults to no retries.
 -startFromScratch,--start_from_scratch                                    Runs all command line functions even if the
                                                                           outputs were previously output successfully.
 -keepIntermediates,--keep_intermediate_outputs                            After a successful run keep the outputs of
                                                                           any Function marked as intermediate.
 -statusTo,--status_email_to &lt;status_email_to&gt;                             Email address to send emails to upon
                                                                           completion or on error.
 -statusFrom,--status_email_from &lt;status_email_from&gt;                       Email address to send emails from upon
                                                                           completion or on error.
 -dot,--dot_graph &lt;dot_graph&gt;                                              Outputs the queue graph to a .dot file.  See:
                                                                           http://en.wikipedia.org/wiki/DOT_language
 -expandedDot,--expanded_dot_graph &lt;expanded_dot_graph&gt;                    Outputs the queue graph of scatter gather to
                                                                           a .dot file.  Otherwise overwrites the
                                                                           dot_graph
 -jobPrefix,--job_name_prefix &lt;job_name_prefix&gt;                            Default name prefix for compute farm jobs.
 -jobProject,--job_project &lt;job_project&gt;                                   Default project for compute farm jobs.
 -jobQueue,--job_queue &lt;job_queue&gt;                                         Default queue for compute farm jobs.
 -jobPriority,--job_priority &lt;job_priority&gt;                                Default priority for jobs.
 -memLimit,--default_memory_limit &lt;default_memory_limit&gt;                   Default memory limit for jobs, in gigabytes.
 -runDir,--run_directory &lt;run_directory&gt;                                   Root directory to run functions from.
 -tempDir,--temp_directory &lt;temp_directory&gt;                                Temp directory to pass to functions.
 -jobSGDir,--job_scatter_gather_directory &lt;job_scatter_gather_directory&gt;   Default directory to place scatter gather
                                                                           output for compute farm jobs.
 -emailHost,--emailSmtpHost &lt;emailSmtpHost&gt;                                Email SMTP host. Defaults to localhost.
 -emailPort,--emailSmtpPort &lt;emailSmtpPort&gt;                                Email SMTP port. Defaults to 465 for ssl,
                                                                           otherwise 25.
 -emailTLS,--emailUseTLS                                                   Email should use TLS. Defaults to false.
 -emailSSL,--emailUseSSL                                                   Email should use SSL. Defaults to false.
 -emailUser,--emailUsername &lt;emailUsername&gt;                                Email SMTP username. Defaults to none.
 -emailPassFile,--emailPasswordFile &lt;emailPasswordFile&gt;                    Email SMTP password file. Defaults to none.
 -emailPass,--emailPassword &lt;emailPassword&gt;                                Email SMTP password. Defaults to none. Not
                                                                           secure! See emailPassFile.
 -l,--logging_level &lt;logging_level&gt;                                        Set the minimum level of logging, i.e.
                                                                           setting INFO get's you INFO up to FATAL,
                                                                           setting ERROR gets you ERROR and FATAL level
                                                                           logging.
 -log,--log_to_file &lt;log_to_file&gt;                                          Set the logging location
 -quiet,--quiet_output_mode                                                Set the logging to quiet mode, no output to
                                                                           stdout
 -debug,--debug_mode                                                       Set the logging file string to include a lot
                                                                           of debugging information (SLOW!)
 -h,--help                                                                 Generate this help message

Arguments for ExampleUnifiedGenotyper:
 -R,--referencefile &lt;referencefile&gt;                          The reference file for the bam files.
 -I,--bamfile &lt;bamfile&gt;                                      Bam file to genotype.
 -L,--intervals &lt;intervals&gt;                                  An optional file with a list of intervals to proccess.
 -filter,--filternames &lt;filternames&gt;                         A optional list of filter names.
 -filterExpression,--filterexpressions &lt;filterexpressions&gt;   An optional list of filter expressions.

##### ERROR ------------------------------------------------------------------------------------------
##### ERROR stack trace
org.broadinstitute.sting.commandline.MissingArgumentException:
Argument with name '--bamfile' (-I) is missing.
Argument with name '--referencefile' (-R) is missing.
        at org.broadinstitute.sting.commandline.ParsingEngine.validate(ParsingEngine.java:192)
        at org.broadinstitute.sting.commandline.ParsingEngine.validate(ParsingEngine.java:172)
        at org.broadinstitute.sting.commandline.CommandLineProgram.start(CommandLineProgram.java:199)
        at org.broadinstitute.sting.queue.QCommandLine$.main(QCommandLine.scala:57)
        at org.broadinstitute.sting.queue.QCommandLine.main(QCommandLine.scala)
##### ERROR ------------------------------------------------------------------------------------------
##### ERROR A GATK RUNTIME ERROR has occurred (version 1.0.5504):
##### ERROR
##### ERROR Please visit the wiki to see if this is a known problem
##### ERROR If not, please post the error, with stack trace, to the GATK forum
##### ERROR Visit our wiki for extensive documentation http://www.broadinstitute.org/gsa/wiki
##### ERROR Visit our forum to view answers to commonly asked questions http://getsatisfaction.com/gsa
##### ERROR
##### ERROR MESSAGE: Argument with name '--bamfile' (-I) is missing.
##### ERROR Argument with name '--referencefile' (-R) is missing.
##### ERROR ------------------------------------------------------------------------------------------</code class="pre_md"></pre>
<p>To dry run the pipeline:</p>
<pre><code class="pre_md">java \
  -Djava.io.tmpdir=tmp \
  -jar dist/Queue.jar \
  -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala \
  -R human_b36_both.fasta \
  -I pilot2_daughters.chr20.10k-11k.bam \
  -L chr20.interval_list \
  -filter StrandBias -filterExpression "SB&gt;=0.10" \
  -filter AlleleBalance -filterExpression "AB&gt;=0.75" \
  -filter QualByDepth -filterExpression "QD&lt;5" \
  -filter HomopolymerRun -filterExpression "HRun&gt;=4"</code class="pre_md"></pre>
<p>The dry run output should appear similar to this:</p>
<pre><code class="pre_md">INFO  10:45:00,354 QScriptManager - Compiling 1 QScript
INFO  10:45:04,855 QScriptManager - Compilation complete
INFO  10:45:05,058 HelpFormatter - ---------------------------------------------------------
INFO  10:45:05,059 HelpFormatter - Program Name: org.broadinstitute.sting.queue.QCommandLine
INFO  10:45:05,059 HelpFormatter - Program Args: -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala -R human_b36_both.fasta -I pilot2_daughters.chr20.10k-11k.bam -L chr20.interval_list -filter StrandBias -filterExpression SB&gt;=0.10 -filter AlleleBalance -filterExpression AB&gt;=0.75 -filter QualByDepth -filterExpression QD&lt;5 -filter HomopolymerRun -filterExpression HRun&gt;=4 
INFO  10:45:05,059 HelpFormatter - Date/Time: 2011/03/24 10:45:05
INFO  10:45:05,059 HelpFormatter - ---------------------------------------------------------
INFO  10:45:05,059 HelpFormatter - ---------------------------------------------------------
INFO  10:45:05,061 QCommandLine - Scripting ExampleUnifiedGenotyper
INFO  10:45:05,150 QCommandLine - Added 4 functions
INFO  10:45:05,150 QGraph - Generating graph.
INFO  10:45:05,169 QGraph - Generating scatter gather jobs.
INFO  10:45:05,182 QGraph - Removing original jobs.
INFO  10:45:05,183 QGraph - Adding scatter gather jobs.
INFO  10:45:05,231 QGraph - Regenerating graph.
INFO  10:45:05,247 QGraph - -------
INFO  10:45:05,252 QGraph - Pending: IntervalScatterFunction /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-1/scatter.intervals /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-2/scatter.intervals /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-3/scatter.intervals
INFO  10:45:05,253 QGraph - Log: /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/scatter/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,254 QGraph - -------
INFO  10:45:05,279 QGraph - Pending: java -Xmx2g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -I /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.bam -L /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-1/scatter.intervals -R /Users/kshakir/src/Sting/human_b36_both.fasta -o /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-1/pilot2_daughters.chr20.10k-11k.unfiltered.vcf
INFO  10:45:05,279 QGraph - Log: /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-1/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,279 QGraph - -------
INFO  10:45:05,283 QGraph - Pending: java -Xmx2g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -I /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.bam -L /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-2/scatter.intervals -R /Users/kshakir/src/Sting/human_b36_both.fasta -o /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-2/pilot2_daughters.chr20.10k-11k.unfiltered.vcf
INFO  10:45:05,283 QGraph - Log: /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-2/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,283 QGraph - -------
INFO  10:45:05,287 QGraph - Pending: java -Xmx2g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -I /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.bam -L /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-3/scatter.intervals -R /Users/kshakir/src/Sting/human_b36_both.fasta -o /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-3/pilot2_daughters.chr20.10k-11k.unfiltered.vcf
INFO  10:45:05,287 QGraph - Log: /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-3/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,288 QGraph - -------
INFO  10:45:05,288 QGraph - Pending: SimpleTextGatherFunction /Users/kshakir/src/Sting/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,288 QGraph - Log: /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/gather-jobOutputFile/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,289 QGraph - -------
INFO  10:45:05,291 QGraph - Pending: java -Xmx1g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T CombineVariants -L /Users/kshakir/src/Sting/chr20.interval_list -R /Users/kshakir/src/Sting/human_b36_both.fasta -B:input0,VCF /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-1/pilot2_daughters.chr20.10k-11k.unfiltered.vcf -B:input1,VCF /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-2/pilot2_daughters.chr20.10k-11k.unfiltered.vcf -B:input2,VCF /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/temp-3/pilot2_daughters.chr20.10k-11k.unfiltered.vcf -o /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.unfiltered.vcf -priority input0,input1,input2 -assumeIdenticalSamples
INFO  10:45:05,291 QGraph - Log: /Users/kshakir/src/Sting/queueScatterGather/Q-60018@bmef8-d8e-1-sg/gather-out/Q-60018@bmef8-d8e-1.out
INFO  10:45:05,292 QGraph - -------
INFO  10:45:05,296 QGraph - Pending: java -Xmx2g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T VariantEval -L /Users/kshakir/src/Sting/chr20.interval_list -R /Users/kshakir/src/Sting/human_b36_both.fasta -B:eval,VCF /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.unfiltered.vcf -o /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.unfiltered.eval
INFO  10:45:05,296 QGraph - Log: /Users/kshakir/src/Sting/Q-60018@bmef8-d8e-2.out
INFO  10:45:05,296 QGraph - -------
INFO  10:45:05,299 QGraph - Pending: java -Xmx2g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T VariantFiltration -L /Users/kshakir/src/Sting/chr20.interval_list -R /Users/kshakir/src/Sting/human_b36_both.fasta -B:vcf,VCF /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.unfiltered.vcf -o /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.filtered.vcf -filter SB&gt;=0.10 -filter AB&gt;=0.75 -filter QD&lt;5 -filter HRun&gt;=4 -filterName StrandBias -filterName AlleleBalance -filterName QualByDepth -filterName HomopolymerRun
INFO  10:45:05,299 QGraph - Log: /Users/kshakir/src/Sting/Q-60018@bmef8-d8e-3.out
INFO  10:45:05,302 QGraph - -------
INFO  10:45:05,303 QGraph - Pending: java -Xmx2g -Djava.io.tmpdir=/Users/kshakir/src/Sting/tmp -cp "/Users/kshakir/src/Sting/dist/Queue.jar" org.broadinstitute.sting.gatk.CommandLineGATK -T VariantEval -L /Users/kshakir/src/Sting/chr20.interval_list -R /Users/kshakir/src/Sting/human_b36_both.fasta -B:eval,VCF /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.filtered.vcf -o /Users/kshakir/src/Sting/pilot2_daughters.chr20.10k-11k.filtered.eval
INFO  10:45:05,303 QGraph - Log: /Users/kshakir/src/Sting/Q-60018@bmef8-d8e-4.out
INFO  10:45:05,304 QGraph - Dry run completed successfully!
INFO  10:45:05,304 QGraph - Re-run with "-run" to execute the functions.
INFO  10:45:05,304 QCommandLine - Done</code class="pre_md"></pre>
<h3>8. Using traits to pass common values between QScripts to CommandLineFunctions</h3>
<p><code>QScript</code> files often create multiple <code>CommandLineFunctions</code> with similar arguments. Use various scala tricks such as inner classes, traits / mixins, etc. to reuse variables.</p>
<ul>
<li>
<p>A <code>self type</code> can be useful to distinguish between <code>this</code>. We use <code>qscript</code> as an alias for the QScript's <code>this</code> to distinguish from the <code>this</code> inside of inner classes or traits.</p>
</li>
<li>A <code>trait mixin</code> can be used to reuse functionality. The trait below is designed to copy values from the QScript and then is mixed into different instances of the functions.</li>
</ul>
<p>See the following example:</p>
<pre><code class="pre_md">class MyScript extends org.broadinstitute.sting.queue.QScript {
  // Create an alias 'qscript' for 'MyScript.this'
  qscript =&gt;

  // This is a script argument
  @Argument(doc="message to display")
  var message: String = _

  // This is a script argument
  @Argument(doc="number of times to display")
  var count: Int = _

  trait ReusableArguments extends MyCommandLineFunction {
    // Whenever a function is created 'with' this trait, it will copy the message.
    this.commandLineMessage = qscript.message
  }

  abstract class MyCommandLineFunction extends CommandLineFunction {
     // This is a per command line argument
     @Argument(doc="message to display")
     var commandLineMessage: String = _
  }

  class MyEchoFunction extends MyCommandLineFunction {
     def commandLine = "echo " + commandLineMessage
  }

  class MyAlsoEchoFunction extends MyCommandLineFunction {
     def commandLine = "echo also " + commandLineMessage
  }

  def script = {
    for (i &lt;- 1 to count) {
      val echo = new MyEchoFunction with ReusableArguments
      val alsoEcho = new MyAlsoEchoFunction with ReusableArguments
      add(echo, alsoEcho)
    }
  }
}</code class="pre_md"></pre>
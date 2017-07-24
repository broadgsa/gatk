## (howto) Run Queue for the first time

http://gatkforums.broadinstitute.org/gatk/discussion/1288/howto-run-queue-for-the-first-time

<h4>Objective</h4>
<p>Run a basic analysis command on example data, parallelized with Queue.</p>
<h4>Prerequisites</h4>
<ul>
<li>Successfully completed <a href="http://gatkforums.broadinstitute.org/discussion/1287/how-to-test-your-queue-installation">&quot;How to test your Queue installation&quot;</a> and <a href="http://gatkforums.broadinstitute.org/discussion/1209/how-to-run-the-gatk-for-the-first-time#latest">&quot;How to run GATK for the first time&quot;</a></li>
<li><a href="http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it">GATK resource bundle</a> downloaded </li>
</ul>
<h4>Steps</h4>
<ol>
<li>Set up a dry run of Queue</li>
<li>Run the analysis for real</li>
<li>Running on a computing farm</li>
</ol>
<hr />
<h3>1. Set up a dry run of Queue</h3>
<p>One very cool feature of Queue is that you can test your script by doing a &quot;dry run&quot;. That means Queue will prepare the analysis and build the scatter commands, but not actually run them. This makes it easier to check the sanity of your script and command. </p>
<p>Here we're going to set up a dry run of a CountReads analysis. You should be familiar with the CountReads walker and the example files from the bundles, as used in the basic &quot;GATK for the first time&quot; tutorial. In addition, we're going to use the example QScript called <code>ExampleCountReads.scala</code> provided in the Queue package download. </p>
<h4>Action</h4>
<p>Type the following command:</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar Queue.jar -S ExampleCountReads.scala -R exampleFASTA.fasta -I exampleBAM.bam</code class="pre_md"></pre>
<p>where <code>-S ExampleCountReads.scala</code> specifies which QScript we want to run, <code>-R exampleFASTA.fasta</code> specifies the reference sequence, and <code>-I exampleBAM.bam</code> specifies the file of aligned reads we want to analyze.</p>
<h4>Expected Result</h4>
<p>After a few seconds you should see output that looks nearly identical to this:</p>
<pre><code class="pre_md">INFO  00:30:45,527 QScriptManager - Compiling 1 QScript 
INFO  00:30:52,869 QScriptManager - Compilation complete 
INFO  00:30:53,284 HelpFormatter - ---------------------------------------------------------------------- 
INFO  00:30:53,284 HelpFormatter - Queue v2.0-36-gf5c1c1a, Compiled 2012/08/08 20:18:21 
INFO  00:30:53,284 HelpFormatter - Copyright (c) 2012 The Broad Institute 
INFO  00:30:53,284 HelpFormatter - Fro support and documentation go to http://www.broadinstitute.org/gatk 
INFO  00:30:53,285 HelpFormatter - Program Args: -S ExampleCountReads.scala -R exampleFASTA.fasta -I exampleBAM.bam 
INFO  00:30:53,285 HelpFormatter - Date/Time: 2012/08/09 00:30:53 
INFO  00:30:53,285 HelpFormatter - ---------------------------------------------------------------------- 
INFO  00:30:53,285 HelpFormatter - ---------------------------------------------------------------------- 
INFO  00:30:53,290 QCommandLine - Scripting ExampleCountReads 
INFO  00:30:53,364 QCommandLine - Added 1 functions 
INFO  00:30:53,364 QGraph - Generating graph. 
INFO  00:30:53,388 QGraph - ------- 
INFO  00:30:53,402 QGraph - Pending:  'java'  '-Xmx1024m'  '-Djava.io.tmpdir=/Users/vdauwera/sandbox/Q2/resources/tmp'  '-cp' '/Users/vdauwera/sandbox/Q2/Queue.jar'  'org.broadinstitute.sting.gatk.CommandLineGATK'  '-T' 'CountReads'  '-I' '/Users/vdauwera/sandbox/Q2/resources/exampleBAM.bam'  '-R' '/Users/vdauwera/sandbox/Q2/resources/exampleFASTA.fasta'  
INFO  00:30:53,403 QGraph - Log:     /Users/vdauwera/sandbox/Q2/resources/ExampleCountReads-1.out 
INFO  00:30:53,403 QGraph - Dry run completed successfully! 
INFO  00:30:53,404 QGraph - Re-run with "-run" to execute the functions. 
INFO  00:30:53,409 QCommandLine - Script completed successfully with 1 total jobs 
INFO  00:30:53,410 QCommandLine - Writing JobLogging GATKReport to file /Users/vdauwera/sandbox/Q2/resources/ExampleCountReads.jobreport.txt </code class="pre_md"></pre>
<p>If you don't see this, check your spelling (GATK commands are case-sensitive), check that the files are in your working directory, and if necessary, re-check that the GATK and Queue are properly installed.</p>
<p>If you do see this output, congratulations! You just successfully ran you first Queue dry run! </p>
<hr />
<h3>2. Run the analysis for real</h3>
<p>Once you have verified that the Queue functions have been generated successfully, you can execute the pipeline by appending <code>-run</code> to the command line.</p>
<h4>Action</h4>
<p>Instead of this command, which we used earlier:</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar Queue.jar -S ExampleCountReads.scala -R exampleFASTA.fasta -I exampleBAM.bam</code class="pre_md"></pre>
<p>this time you type this:</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar Queue.jar -S ExampleCountReads.scala -R exampleFASTA.fasta -I exampleBAM.bam -run</code class="pre_md"></pre>
<p>See the difference?</p>
<h4>Result</h4>
<p>You should see output that looks nearly identical to this:</p>
<pre><code class="pre_md">INFO  00:56:33,688 QScriptManager - Compiling 1 QScript 
INFO  00:56:39,327 QScriptManager - Compilation complete 
INFO  00:56:39,487 HelpFormatter - ---------------------------------------------------------------------- 
INFO  00:56:39,487 HelpFormatter - Queue v2.0-36-gf5c1c1a, Compiled 2012/08/08 20:18:21 
INFO  00:56:39,488 HelpFormatter - Copyright (c) 2012 The Broad Institute 
INFO  00:56:39,488 HelpFormatter - Fro support and documentation go to http://www.broadinstitute.org/gatk 
INFO  00:56:39,489 HelpFormatter - Program Args: -S ExampleCountReads.scala -R exampleFASTA.fasta -I exampleBAM.bam -run 
INFO  00:56:39,490 HelpFormatter - Date/Time: 2012/08/09 00:56:39 
INFO  00:56:39,490 HelpFormatter - ---------------------------------------------------------------------- 
INFO  00:56:39,491 HelpFormatter - ---------------------------------------------------------------------- 
INFO  00:56:39,498 QCommandLine - Scripting ExampleCountReads 
INFO  00:56:39,569 QCommandLine - Added 1 functions 
INFO  00:56:39,569 QGraph - Generating graph. 
INFO  00:56:39,589 QGraph - Running jobs. 
INFO  00:56:39,623 FunctionEdge - Starting:  'java'  '-Xmx1024m'  '-Djava.io.tmpdir=/Users/vdauwera/sandbox/Q2/resources/tmp'  '-cp' '/Users/vdauwera/sandbox/Q2/Queue.jar'  'org.broadinstitute.sting.gatk.CommandLineGATK'  '-T' 'CountReads'  '-I' '/Users/vdauwera/sandbox/Q2/resources/exampleBAM.bam'  '-R' '/Users/vdauwera/sandbox/Q2/resources/exampleFASTA.fasta'  
INFO  00:56:39,623 FunctionEdge - Output written to /Users/GG/codespace/GATK/Q2/resources/ExampleCountReads-1.out 
INFO  00:56:50,301 QGraph - 0 Pend, 1 Run, 0 Fail, 0 Done 
INFO  00:57:09,827 FunctionEdge - Done:  'java'  '-Xmx1024m'  '-Djava.io.tmpdir=/Users/vdauwera/sandbox/Q2/resources/tmp'  '-cp' '/Users/vdauwera/sandbox/Q2/resources/Queue.jar'  'org.broadinstitute.sting.gatk.CommandLineGATK'  '-T' 'CountReads'  '-I' '/Users/vdauwera/sandbox/Q2/resources/exampleBAM.bam'  '-R' '/Users/vdauwera/sandbox/Q2/resources/exampleFASTA.fasta'  
INFO  00:57:09,828 QGraph - 0 Pend, 0 Run, 0 Fail, 1 Done 
INFO  00:57:09,835 QCommandLine - Script completed successfully with 1 total jobs 
INFO  00:57:09,835 QCommandLine - Writing JobLogging GATKReport to file /Users/vdauwera/sandbox/Q2/resources/ExampleCountReads.jobreport.txt 
INFO  00:57:10,107 QCommandLine - Plotting JobLogging GATKReport to file /Users/vdauwera/sandbox/Q2/resources/ExampleCountReads.jobreport.pdf 
WARN  00:57:18,597 RScriptExecutor - RScript exited with 1. Run with -l DEBUG for more info. </code class="pre_md"></pre>
<p>Great! It works!</p>
<p>The results of the traversal will be written to a file in the current directory. The name of the file will be printed in the output, ExampleCountReads.out in this example.</p>
<p>If for some reason the run was interrupted, in most cases you can resume by just launching the command. Queue will pick up where it left off without redoing the parts that ran successfully. </p>
<hr />
<h3>3. Running on a computing farm</h3>
<p>Run with <code>-bsub</code> to run on LSF, or for early Grid Engine support see <a href="http://www.broadinstitute.org/gatk/guide/article?id=1313">Queue with Grid Engine</a>.</p>
<p>See also <a href="http://www.broadinstitute.org/gatk/guide/article?id=1311">QFunction and Command Line Options</a> for more info on Queue options.</p>
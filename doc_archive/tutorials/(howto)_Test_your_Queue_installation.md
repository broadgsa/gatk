## (howto) Test your Queue installation

http://gatkforums.broadinstitute.org/gatk/discussion/1287/howto-test-your-queue-installation

<h4>Objective</h4>
<p>Test that Queue is correctly installed, and that the supporting tools like Java are in your path.</p>
<h4>Prerequisites</h4>
<ul>
<li>Basic familiarity with the command-line environment</li>
<li>Understand what is a PATH variable</li>
<li>GATK installed</li>
<li>Queue downloaded and placed on path</li>
</ul>
<h4>Steps</h4>
<ol>
<li>Invoke the Queue usage/help message</li>
<li>Troubleshooting</li>
</ol>
<hr />
<h3>1. Invoke the Queue usage/help message</h3>
<p>The command we're going to run is a very simple command that asks Queue to print out a list of available command-line arguments and options. It is so simple that it will ALWAYS work if your Queue package is installed correctly.</p>
<p>Note that this command is also helpful when you're trying to remember something like the right spelling or short name for an argument and for whatever reason you don't have access to the web-based documentation.  </p>
<h4>Action</h4>
<p>Type the following command:</p>
<pre><code class="pre_md">java -jar &lt;path to Queue.jar&gt; --help</code class="pre_md"></pre>
<p>replacing the <code>&lt;path to Queue.jar&gt;</code> bit with the path you have set up in your command-line environment.</p>
<h4>Expected Result</h4>
<p>You should see usage output similar to the following:</p>
<pre><code class="pre_md">usage: java -jar Queue.jar -S &lt;script&gt; [-jobPrefix &lt;job_name_prefix&gt;] [-jobQueue &lt;job_queue&gt;] [-jobProject &lt;job_project&gt;]
       [-jobSGDir &lt;job_scatter_gather_directory&gt;] [-memLimit &lt;default_memory_limit&gt;] [-runDir &lt;run_directory&gt;] [-tempDir
       &lt;temp_directory&gt;] [-emailHost &lt;emailSmtpHost&gt;] [-emailPort &lt;emailSmtpPort&gt;] [-emailTLS] [-emailSSL] [-emailUser
       &lt;emailUsername&gt;] [-emailPass &lt;emailPassword&gt;] [-emailPassFile &lt;emailPasswordFile&gt;] [-bsub] [-run] [-dot &lt;dot_graph&gt;]
       [-expandedDot &lt;expanded_dot_graph&gt;] [-startFromScratch] [-status] [-statusFrom &lt;status_email_from&gt;] [-statusTo
       &lt;status_email_to&gt;] [-keepIntermediates] [-retry &lt;retry_failed&gt;] [-l &lt;logging_level&gt;] [-log &lt;log_to_file&gt;] [-quiet]
       [-debug] [-h]

 -S,--script &lt;script&gt;                                                      QScript scala file
 -jobPrefix,--job_name_prefix &lt;job_name_prefix&gt;                            Default name prefix for compute farm jobs.
 -jobQueue,--job_queue &lt;job_queue&gt;                                         Default queue for compute farm jobs.
 -jobProject,--job_project &lt;job_project&gt;                                   Default project for compute farm jobs.
 -jobSGDir,--job_scatter_gather_directory &lt;job_scatter_gather_directory&gt;   Default directory to place scatter gather
                                                                           output for compute farm jobs.
 -memLimit,--default_memory_limit &lt;default_memory_limit&gt;                   Default memory limit for jobs, in gigabytes.
 -runDir,--run_directory &lt;run_directory&gt;                                   Root directory to run functions from.
 -tempDir,--temp_directory &lt;temp_directory&gt;                                Temp directory to pass to functions.
 -emailHost,--emailSmtpHost &lt;emailSmtpHost&gt;                                Email SMTP host. Defaults to localhost.
 -emailPort,--emailSmtpPort &lt;emailSmtpPort&gt;                                Email SMTP port. Defaults to 465 for ssl,
                                                                           otherwise 25.
 -emailTLS,--emailUseTLS                                                   Email should use TLS. Defaults to false.
 -emailSSL,--emailUseSSL                                                   Email should use SSL. Defaults to false.
 -emailUser,--emailUsername &lt;emailUsername&gt;                                Email SMTP username. Defaults to none.
 -emailPass,--emailPassword &lt;emailPassword&gt;                                Email SMTP password. Defaults to none. Not
                                                                           secure! See emailPassFile.
 -emailPassFile,--emailPasswordFile &lt;emailPasswordFile&gt;                    Email SMTP password file. Defaults to none.
 -bsub,--bsub_all_jobs                                                     Use bsub to submit jobs
 -run,--run_scripts                                                        Run QScripts.  Without this flag set only
                                                                           performs a dry run.
 -dot,--dot_graph &lt;dot_graph&gt;                                              Outputs the queue graph to a .dot file.  See:
                                                                           http://en.wikipedia.org/wiki/DOT_language
 -expandedDot,--expanded_dot_graph &lt;expanded_dot_graph&gt;                    Outputs the queue graph of scatter gather to
                                                                           a .dot file.  Otherwise overwrites the
                                                                           dot_graph
 -startFromScratch,--start_from_scratch                                    Runs all command line functions even if the
                                                                           outputs were previously output successfully.
 -status,--status                                                          Get status of jobs for the qscript
 -statusFrom,--status_email_from &lt;status_email_from&gt;                       Email address to send emails from upon
                                                                           completion or on error.
 -statusTo,--status_email_to &lt;status_email_to&gt;                             Email address to send emails to upon
                                                                           completion or on error.
 -keepIntermediates,--keep_intermediate_outputs                            After a successful run keep the outputs of
                                                                           any Function marked as intermediate.
 -retry,--retry_failed &lt;retry_failed&gt;                                      Retry the specified number of times after a
                                                                           command fails.  Defaults to no retries.
 -l,--logging_level &lt;logging_level&gt;                                        Set the minimum level of logging, i.e.
                                                                           setting INFO get's you INFO up to FATAL,
                                                                           setting ERROR gets you ERROR and FATAL level
                                                                           logging.
 -log,--log_to_file &lt;log_to_file&gt;                                          Set the logging location
 -quiet,--quiet_output_mode                                                Set the logging to quiet mode, no output to
                                                                           stdout
 -debug,--debug_mode                                                       Set the logging file string to include a lot
                                                                           of debugging information (SLOW!)
 -h,--help                                                                 Generate this help message</code class="pre_md"></pre>
<p>If you see this message, your Queue installation is ok. You're good to go! If you don't see this message, and instead get an error message, proceed to the next section on troubleshooting.  </p>
<hr />
<h3>2. Troubleshooting</h3>
<p>Let's try to figure out what's not working.  </p>
<h4>Action</h4>
<p>First, make sure that your Java version is at least 1.6, by typing the following command:</p>
<pre><code class="pre_md">java -version</code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>You should see something similar to the following text:</p>
<pre><code class="pre_md">java version "1.6.0_12"
Java(TM) SE Runtime Environment (build 1.6.0_12-b04)
Java HotSpot(TM) 64-Bit Server VM (build 11.2-b01, mixed mode)  </code class="pre_md"></pre>
<h4>Remedial actions</h4>
<p>If the version is less then 1.6, install the newest version of Java onto the system. If you instead see something like </p>
<pre><code class="pre_md">java: Command not found  </code class="pre_md"></pre>
<p>make sure that java is installed on your machine, and that your PATH variable contains the path to the java executables. </p>
<p>On a Mac running OS X 10.5+, you may need to run /Applications/Utilities/Java Preferences.app and drag Java SE 6 to the top to make your machine run version 1.6, even if it has been installed.</p>
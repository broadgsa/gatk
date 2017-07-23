## QFunction and Command Line Options

http://gatkforums.broadinstitute.org/gatk/discussion/1311/qfunction-and-command-line-options

<p>These are the most popular Queue command line options. For a complete and up to date list run with <code>--help</code> or <code>-h</code>. QScripts may also add additional command line options.</p>
<p><strong>Please note that this page is out of date. We hope to update it in future but have no resources to do so at present. If you run into trouble using any of the command line arguments listed here, we recommend you check the source code for the Q arguments <a href="https://github.com/broadgsa/gatk/blob/master/public/gatk-queue/src/main/scala/org/broadinstitute/gatk/queue/QSettings.scala">here</a>. Apologies for the inconvenience.</strong></p>
<hr />
<h3>1. Queue Command Line Options</h3>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Command Line Argument</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Default</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>-run</code></td>
<td style="text-align: left;">If passed the scripts are run. If not passed a dry run is executed.</td>
<td style="text-align: left;">dry run</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobRunner &lt;jobrunner&gt;</code></td>
<td style="text-align: left;">The job runner to dispatch jobs. Setting to <code>Lsf706</code>, <code>GridEngine</code>, or <code>Drmaa</code> will dispatch jobs to LSF or Grid Engine using the job settings (see below). Defaults to <code>Shell</code> which runs jobs on a local shell one at a time.</td>
<td style="text-align: left;"><code>Shell</code></td>
</tr>
<tr>
<td style="text-align: left;"><code>-bsub</code></td>
<td style="text-align: left;">Alias for <code>-jobRunner Lsf706</code></td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-qsub</code></td>
<td style="text-align: left;">Alias for <code>-jobRunner GridEngine</code></td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-status</code></td>
<td style="text-align: left;">Prints out a summary progress. If a QScript is currently running via <code>-run</code>, you can run the same command line with <code>-status</code> instead to print a summary of progress.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-retry &lt;count&gt;</code></td>
<td style="text-align: left;">Retries a QFunction that returns a non-zero exit code up to count times. The QFunction must not have set <code>jobRestartable</code> to <code>false</code>.</td>
<td style="text-align: left;"><code>0</code> = no retries</td>
</tr>
<tr>
<td style="text-align: left;"><code>-startFromScratch</code></td>
<td style="text-align: left;">Restarts the graph from the beginning. If not specified for each output file specified on a QFunction, ex: <code>/path/to/output.file</code>, Queue will not re-run the job if a <code>.done</code> file is found for the all the outputs, ex: <code>/path/to/.output.file.done</code>.</td>
<td style="text-align: left;">use <code>.done</code> files to determine if jobs are complete</td>
</tr>
<tr>
<td style="text-align: left;"><code>-keepIntermediates</code></td>
<td style="text-align: left;">By default Queue deletes the output files of QFunctions that set <code>.isIntermediate</code> to <code>true</code>.</td>
<td style="text-align: left;">delete intermediate files</td>
</tr>
<tr>
<td style="text-align: left;"><code>-statusTo &lt;email&gt;</code></td>
<td style="text-align: left;">Email address to send status to whenever a) A job fails, or b) Queue has run all the functions it can run and is exiting.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-statusFrom &lt;email&gt;</code></td>
<td style="text-align: left;">Email address to send status emails from.</td>
<td style="text-align: left;"><code>user@local.domain</code></td>
</tr>
<tr>
<td style="text-align: left;"><code>-dot &lt;file&gt;</code></td>
<td style="text-align: left;">If set renders the job graph to a <a href="http://en.wikipedia.org/wiki/DOT_language">dot file</a>.</td>
<td style="text-align: left;">not rendered</td>
</tr>
<tr>
<td style="text-align: left;"><code>-l &lt;logging_level&gt;</code></td>
<td style="text-align: left;">The minimum level of logging, <code>DEBUG</code>, <code>INFO</code>, <code>WARN</code>, or <code>FATAL</code>.</td>
<td style="text-align: left;"><code>INFO</code></td>
</tr>
<tr>
<td style="text-align: left;"><code>-log &lt;file&gt;</code></td>
<td style="text-align: left;">Sets the location to save log output in addition to standard out.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-debug</code></td>
<td style="text-align: left;">Set the logging to include a lot of debugging information (SLOW!)</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobReport</code></td>
<td style="text-align: left;">Path to write the job report text file. If R is installed and available on the <code>$PATH</code> then a pdf will be generated visualizing the job report.</td>
<td style="text-align: left;"><code>jobPrefix.jobreport.txt</code></td>
</tr>
<tr>
<td style="text-align: left;"><code>-disableJobReport</code></td>
<td style="text-align: left;">Disables writing the job report.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-help</code></td>
<td style="text-align: left;">Lists all of the command line arguments with their descriptions.</td>
<td style="text-align: left;">not set</td>
</tr>
</tbody>
</table>
<h3>2. QFunction Options</h3>
<p>The following options can be specified on the command line over overridden per QFunction.</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Command Line Argument</th>
<th style="text-align: left;">QFunction Property</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Default</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>-jobPrefix</code></td>
<td style="text-align: left;"><code>.jobName</code></td>
<td style="text-align: left;">The unique name of the job. Used to prefix directories and log files. Use <code>-jobNamePrefix</code> on the Queue command line to replace the default prefix <code>Q-&lt;processid&gt;@&lt;host&gt;</code>.</td>
<td style="text-align: left;"><code>&lt;jobNamePrefix&gt;-&lt;jobNumber&gt;</code></td>
</tr>
<tr>
<td style="text-align: left;">N/A</td>
<td style="text-align: left;"><code>.jobOutputFile</code></td>
<td style="text-align: left;">Captures <code>stdout</code> and if <code>jobErrorFile</code> is null it captures <code>stderr</code> as well.</td>
<td style="text-align: left;"><code>&lt;jobName&gt;.out</code></td>
</tr>
<tr>
<td style="text-align: left;">N/A</td>
<td style="text-align: left;"><code>.jobErrorFile</code></td>
<td style="text-align: left;">If not null captures stderr.</td>
<td style="text-align: left;"><code>null</code></td>
</tr>
<tr>
<td style="text-align: left;">N/A</td>
<td style="text-align: left;"><code>.commandDirectory</code></td>
<td style="text-align: left;">The directory to execute the command line from.</td>
<td style="text-align: left;">current directory</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobProject</code></td>
<td style="text-align: left;"><code>.jobProject</code></td>
<td style="text-align: left;">The project name for the job.</td>
<td style="text-align: left;">default job runner project</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobQueue</code></td>
<td style="text-align: left;"><code>.jobQueue</code></td>
<td style="text-align: left;">The queue to dispatch the job.</td>
<td style="text-align: left;">default job runner queue</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobPriority</code></td>
<td style="text-align: left;"><code>.jobPriority</code></td>
<td style="text-align: left;">The dispatch priority for the job. Lowest priority = <code>0</code>. Highest priority = <code>100</code>.</td>
<td style="text-align: left;">default job runner priority</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobNative</code></td>
<td style="text-align: left;"><code>.jobNativeArgs</code></td>
<td style="text-align: left;">Native args to pass to the job runner. Currently only supported in GridEngine and Drmaa. The string is concatenated to the native arguments passed over DRMAA. Example: <code>-w n</code>.</td>
<td style="text-align: left;">none</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobResReq</code></td>
<td style="text-align: left;"><code>.jobResourceRequests</code></td>
<td style="text-align: left;">Resource requests to pass to the job runner. On GridEngine this is multiple <code>-l &lt;req&gt;</code>. On LSF a single <code>-R &lt;req&gt;</code> is generated.</td>
<td style="text-align: left;">memory reservations and limits on LSF and GridEngine</td>
</tr>
<tr>
<td style="text-align: left;"><code>-jobEnv</code></td>
<td style="text-align: left;"><code>.jobEnvironmentNames</code></td>
<td style="text-align: left;">Predefined environment names to pass to the job runner. On GridEngine this is <code>-pe &lt;env&gt;</code>. On LSF this is <code>-a &lt;env&gt;</code>.</td>
<td style="text-align: left;">none</td>
</tr>
<tr>
<td style="text-align: left;"><code>-memLimit</code></td>
<td style="text-align: left;"><code>.memoryLimit</code></td>
<td style="text-align: left;">The memory limit for the job in gigabytes. Used to populate the variables residentLimit and residentRequest which can also be set separately.</td>
<td style="text-align: left;">default job runner memory limit</td>
</tr>
<tr>
<td style="text-align: left;"><code>-resMemLimit</code></td>
<td style="text-align: left;"><code>.residentLimit</code></td>
<td style="text-align: left;">Limit for the resident memory in gigabytes. On GridEngine this is <code>-l mem_free=&lt;mem&gt;</code>. On LSF this is <code>-R rusage[mem=&lt;mem&gt;]</code>.</td>
<td style="text-align: left;"><code>memoryLimit</code> * 1.2</td>
</tr>
<tr>
<td style="text-align: left;"><code>-resMemReq</code></td>
<td style="text-align: left;"><code>.residentRequest</code></td>
<td style="text-align: left;">Requested amount of resident memory in gigabytes. On GridEngine this is <code>-l h_rss=&lt;mem&gt;</code>. On LSF this is <code>-R rusage[select=&lt;mem&gt;]</code>.</td>
<td style="text-align: left;"><code>memoryLimit</code></td>
</tr>
</tbody>
</table>
<h3>3. Email Status Options</h3>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Command Line Argument</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Default</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>-emailHost &lt;hostname&gt;</code></td>
<td style="text-align: left;">SMTP host name</td>
<td style="text-align: left;">localhost</td>
</tr>
<tr>
<td style="text-align: left;"><code>-emailPort &lt;port&gt;</code></td>
<td style="text-align: left;">SMTP port</td>
<td style="text-align: left;">25</td>
</tr>
<tr>
<td style="text-align: left;"><code>-emailTLS</code></td>
<td style="text-align: left;">If set uses TLS.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-emailSSL</code></td>
<td style="text-align: left;">If set uses SSL.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-emailUser &lt;username&gt;</code></td>
<td style="text-align: left;">If set along with emailPass or emailPassFile authenticates the email with this username.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-emailPassFile &lt;file&gt;</code></td>
<td style="text-align: left;">If emailUser is also set authenticates the email with contents of the file.</td>
<td style="text-align: left;">not set</td>
</tr>
<tr>
<td style="text-align: left;"><code>-emailPass &lt;password&gt;</code></td>
<td style="text-align: left;">If emailUser is also set authenticates the email with this password. NOT SECURE: Use emailPassFile instead!</td>
<td style="text-align: left;">not set</td>
</tr>
</tbody>
</table>
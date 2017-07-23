## Queue custom job schedulers

http://gatkforums.broadinstitute.org/gatk/discussion/1347/queue-custom-job-schedulers

<h2>Implementing a Queue JobRunner</h2>
<p>The following scala methods need to be implemented for a new JobRunner. See the implementations of <a href="https://github.com/broadgsa/gatk/blob/master/public/queue-framework/src/main/scala/org/broadinstitute/sting/queue/engine/gridengine/GridEngineJobRunner.scala">GridEngine</a> and <a href="https://github.com/broadgsa/gatk/blob/master/public/queue-framework/src/main/scala/org/broadinstitute/sting/queue/engine/lsf/Lsf706JobRunner.scala">LSF</a> for concrete full examples. </p>
<h3>1. class JobRunner.start()</h3>
<p>Start should to copy the settings from the CommandLineFunction into your job scheduler and invoke the command via <code>sh &lt;jobScript&gt;</code>. As an example of what needs to be implemented, here is the current contents of the <code>start()</code> method in <code>MyCustomJobRunner</code> which contains the pseudo code.</p>
<pre><code class="pre_md">  def start() {
    // TODO: Copy settings from function to your job scheduler syntax.

    val mySchedulerJob = new ...

    // Set the display name to 4000 characters of the description (or whatever your max is)
    mySchedulerJob.displayName = function.description.take(4000)

    // Set the output file for stdout
    mySchedulerJob.outputFile = function.jobOutputFile.getPath

    // Set the current working directory
    mySchedulerJob.workingDirectory = function.commandDirectory.getPath

    // If the error file is set specify the separate output for stderr
    if (function.jobErrorFile != null) {
      mySchedulerJob.errFile = function.jobErrorFile.getPath
    }

    // If a project name is set specify the project name
    if (function.jobProject != null) {
      mySchedulerJob.projectName = function.jobProject
    }

    // If the job queue is set specify the job queue
    if (function.jobQueue != null) {
      mySchedulerJob.queue = function.jobQueue
    }

    // If the resident set size is requested pass on the memory request
    if (residentRequestMB.isDefined) {
      mySchedulerJob.jobMemoryRequest = "%dM".format(residentRequestMB.get.ceil.toInt)
    }

    // If the resident set size limit is defined specify the memory limit
    if (residentLimitMB.isDefined) {
      mySchedulerJob.jobMemoryLimit = "%dM".format(residentLimitMB.get.ceil.toInt)
    }

    // If the priority is set (user specified Int) specify the priority
    if (function.jobPriority.isDefined) {
      mySchedulerJob.jobPriority = function.jobPriority.get
    }

    // Instead of running the function.commandLine, run "sh &lt;jobScript&gt;"
    mySchedulerJob.command = "sh " + jobScript

    // Store the status so it can be returned in the status method.
    myStatus = RunnerStatus.RUNNING

    // Start the job and store the id so it can be killed in tryStop
    myJobId = mySchedulerJob.start()
  }</code class="pre_md"></pre>
<h3>2. class JobRunner.status</h3>
<p>The status method should return one of the enum values from <code>org.broadinstitute.sting.queue.engine.RunnerStatus</code>:</p>
<ul>
<li><code>RunnerStatus.RUNNING</code></li>
<li><code>RunnerStatus.DONE</code></li>
<li><code>RunnerStatus.FAILED</code></li>
</ul>
<h3>3. object JobRunner.init()</h3>
<p>Add any initialization code to the companion object static initializer. See the LSF or GridEngine implementations for how this is done.</p>
<h3>4. object JobRunner.tryStop()</h3>
<p>The jobs that are still in <code>RunnerStatus.RUNNING</code> will be passed into this function. <code>tryStop()</code> should send these jobs the equivalent of a <code>Ctrl-C</code> or <code>SIGTERM(15)</code>, or worst case a <code>SIGKILL(9)</code> if <code>SIGTERM</code> is not available.</p>
<h2>Running Queue with a new JobRunner</h2>
<p>Once there is a basic implementation, you can try out the Hello World example with <code>-jobRunner MyJobRunner</code>.</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar dist/Queue.jar -S scala/qscript/examples/HelloWorld.scala -jobRunner MyJobRunner -run</code class="pre_md"></pre>
<p>If all goes well Queue should dispatch the job to your job scheduler and wait until the status returns <code>RunningStatus.DONE</code> and <code>hello world</code> should be echo'ed into the output file, possibly with other log messages.</p>
<p>See [QFunction and Command Line Options]() for more info on Queue options.</p>
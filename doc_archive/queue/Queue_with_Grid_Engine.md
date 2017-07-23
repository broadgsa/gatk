## Queue with Grid Engine

http://gatkforums.broadinstitute.org/gatk/discussion/1313/queue-with-grid-engine

<h3>1. Background</h3>
<p>Thanks to contributions from the community, Queue contains a job runner compatible with Grid Engine 6.2u5.</p>
<p>As of July 2011 this is the currently known list of forked distributions of Sun's Grid Engine 6.2u5. As long as they are <a href="http://gridscheduler.sourceforge.net/javadocs/">JDRMAA 1.0 source compatible</a> with Grid Engine 6.2u5, the compiled Queue code should run against each of these distributions. However we have yet to receive confirmation that Queue works on any of these setups.</p>
<ul>
<li><a href="http://wikis.sun.com/display/gridengine62u7/Home">Oracle Grid Engine 6.2u7</a></li>
<li><a href="http://gridengine.org">Univa Grid Engine Core 8.0.0</a> </li>
<li><a href="http://www.univa.com/products/grid-engine">Univa Grid Engine  8.0.0</a> </li>
<li><a href="https://arc.liv.ac.uk/SGE">Son of Grid Engine 8.0.0a</a></li>
<li><a href="http://www.rocksclusters.org/">Rocks 5.4</a> (includes a Roll for <a href="http://www.rocksclusters.org/roll-documentation/base/5.4/x8106.html#AEN8149">&quot;SGE V62u5&quot;</a>)</li>
<li><a href="http://gridscheduler.sourceforge.net/">Open Grid Scheduler 6.2u5p2</a></li>
</ul>
<p>Our internal QScript integration tests run the same tests on both LSF 7.0.6 and a Grid Engine 6.2u5 cluster setup on older software released by Sun.</p>
<p>If you run into trouble, please let us know. If you would like to contribute additions or bug fixes please create a fork in our <a href="https://github.com/broadgsa/gatk">github repo</a> where we can review and pull in the patch.</p>
<h2>2. Running Queue with GridEngine</h2>
<p>Try out the Hello World example with <code>-jobRunner GridEngine</code>.</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar dist/Queue.jar -S public/scala/qscript/examples/HelloWorld.scala -jobRunner GridEngine -run</code class="pre_md"></pre>
<p>If all goes well Queue should dispatch the job to Grid Engine and wait until the status returns <code>RunningStatus.DONE</code> and &quot;<code>hello world</code> should be echoed into the output file, possibly with other grid engine log messages.</p>
<p>See <a href="http://gatkforums.broadinstitute.org/discussion/1311/qfunction-and-command-line-options">QFunction and Command Line Options</a> for more info on Queue options.</p>
<h2>3. Debugging issues with Queue and GridEngine</h2>
<p>If you run into an error with Queue submitting jobs to GridEngine, first try submitting the HelloWorld example with <code>-memLimit 2</code>:</p>
<pre><code class="pre_md">java -Djava.io.tmpdir=tmp -jar dist/Queue.jar -S public/scala/qscript/examples/HelloWorld.scala -jobRunner GridEngine -run -memLimit 2</code class="pre_md"></pre>
<p>Then try the following GridEngine qsub commands. They are based on what Queue submits via the API when running the <code>HelloWorld.scala</code> example with and without memory reservations and limits: </p>
<pre><code class="pre_md">qsub -w e -V -b y -N echo_hello_world \
  -o test.out -wd $PWD -j y echo hello world

qsub -w e -V -b y -N echo_hello_world \
  -o test.out -wd $PWD -j y \
  -l mem_free=2048M -l h_rss=2458M echo hello world</code class="pre_md"></pre>
<p>One other thing to check is if there is a memory limit on your cluster. For example try submitting jobs with up to 16G.</p>
<pre><code class="pre_md">qsub -w e -V -b y -N echo_hello_world \
  -o test.out -wd $PWD -j y \
  -l mem_free=4096M -l h_rss=4915M echo hello world

qsub -w e -V -b y -N echo_hello_world \
  -o test.out -wd $PWD -j y \
  -l mem_free=8192M -l h_rss=9830M echo hello world

qsub -w e -V -b y -N echo_hello_world \
  -o test.out -wd $PWD -j y \
  -l mem_free=16384M -l h_rss=19960M echo hello world</code class="pre_md"></pre>
<p>If the above tests pass and GridEngine will still not dispatch jobs submitted by Queue please report the issue to our <a href="http://gatkforums.broadinstitute.org/">support forum</a>.</p>
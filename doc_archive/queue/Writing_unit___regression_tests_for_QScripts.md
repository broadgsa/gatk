## Writing unit / regression tests for QScripts

http://gatkforums.broadinstitute.org/gatk/discussion/1353/writing-unit-regression-tests-for-qscripts

<p>In addition to testing walkers individually, you may want to also run integration tests for your QScript pipelines.</p>
<h2>1. Brief comparison to the Walker integration tests</h2>
<ul>
<li>Pipeline tests should use the standard location for testing data.</li>
<li>Pipeline tests use the same test dependencies.</li>
<li>Pipeline tests which generate MD5 results will have the results stored in the MD5 database]. </li>
<li>Pipeline tests, like QScripts, are written in Scala.</li>
<li>Pipeline tests dry-run under the <code>ant</code> target <code>pipelinetest</code> and run under <code>pipelinetestrun</code>.</li>
<li>Pipeline tests class names must end in <code>PipelineTest</code> to run under the <code>ant</code> target.</li>
<li>Pipeline tests should instantiate a <code>PipelineTestSpec</code> and then run it via <code>PipelineTest.exec()</code>.</li>
</ul>
<h2>2. PipelineTestSpec</h2>
<p>When building up a pipeline test spec specify the following variables for your test.</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;">Variable</th>
<th style="text-align: left;">Type</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>args</code></td>
<td style="text-align: left;"><code>String</code></td>
<td style="text-align: left;">The arguments to pass to the Queue test, ex: <code>-S scala/qscript/examples/HelloWorld.scala</code></td>
</tr>
<tr>
<td style="text-align: left;"><code>jobQueue</code></td>
<td style="text-align: left;"><code>String</code></td>
<td style="text-align: left;">Job Queue to run the test.  Default is <code>null</code> which means use <code>hour</code>.</td>
</tr>
<tr>
<td style="text-align: left;"><code>fileMD5s</code></td>
<td style="text-align: left;"><code>Map[Path, MD5]</code></td>
<td style="text-align: left;">Expected MD5 results for each file path.</td>
</tr>
<tr>
<td style="text-align: left;"><code>expectedException</code></td>
<td style="text-align: left;"><code>classOf[Exception]</code></td>
<td style="text-align: left;">Expected exception from the test.</td>
</tr>
</tbody>
</table>
<h2>3. Example PipelineTest</h2>
<p>The following example runs the <code>ExampleCountLoci</code> QScript on a small bam and verifies that the MD5 result is as expected.</p>
<p>It is checked into the Sting repository under <code>scala/test/org/broadinstitute/sting/queue/pipeline/examples/ExampleCountLociPipelineTest.scala</code></p>
<pre><code class="pre_md">package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}
import org.broadinstitute.sting.BaseTest

class ExampleCountLociPipelineTest {
  @Test
  def testCountLoci {
    val testOut = "count.out"
    val spec = new PipelineTestSpec
    spec.name = "countloci"
    spec.args = Array(
      " -S scala/qscript/examples/ExampleCountLoci.scala",
      " -R " + BaseTest.hg18Reference,
      " -I " + BaseTest.validationDataLocation + "small_bam_for_countloci.bam",
      " -o " + testOut).mkString
    spec.fileMD5s += testOut -&gt; "67823e4722495eb10a5e4c42c267b3a6"
    PipelineTest.executeTest(spec)
  }
}</code class="pre_md"></pre>
<h2>3. Running Pipeline Tests</h2>
<h3>Dry Run</h3>
<p>To test if the script is at least compiling with your arguments run <code>ant pipelinetest</code> specifying the name of your class to <code>-Dsingle</code>:</p>
<pre><code class="pre_md">ant pipelinetest -Dsingle=ExampleCountLociPipelineTest</code class="pre_md"></pre>
<p>Sample output:</p>
<pre><code class="pre_md">   [testng] --------------------------------------------------------------------------------
   [testng] Executing test countloci with Queue arguments: -S scala/qscript/examples/ExampleCountLoci.scala -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/small_bam_for_countloci.bam -o count.out -bsub -l WARN -tempDir pipelinetests/countloci/temp/ -runDir pipelinetests/countloci/run/ -jobQueue hour
   [testng]   =&gt; countloci PASSED DRY RUN
   [testng] PASSED: testCountLoci</code class="pre_md"></pre>
<h3>Run</h3>
<p>As of July 2011 the pipeline tests run against LSF 7.0.6 and Grid Engine 6.2u5. To include these two packages in your environment use the hidden dotkit <code>.combined_LSF_SGE</code>.</p>
<pre><code class="pre_md">reuse .combined_LSF_SGE</code class="pre_md"></pre>
<p>Once you are satisfied that the dry run has completed without error, to actually run the pipeline test run <code>ant pipelinetestrun</code>.</p>
<pre><code class="pre_md">ant pipelinetestrun -Dsingle=ExampleCountLociPipelineTest</code class="pre_md"></pre>
<p>Sample output:</p>
<pre><code class="pre_md">   [testng] --------------------------------------------------------------------------------
   [testng] Executing test countloci with Queue arguments: -S scala/qscript/examples/ExampleCountLoci.scala -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/small_bam_for_countloci.bam -o count.out -bsub -l WARN -tempDir pipelinetests/countloci/temp/ -runDir pipelinetests/countloci/run/ -jobQueue hour -run
   [testng] ##### MD5 file is up to date: integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] Checking MD5 for pipelinetests/countloci/run/count.out [calculated=67823e4722495eb10a5e4c42c267b3a6, expected=67823e4722495eb10a5e4c42c267b3a6]
   [testng]   =&gt; countloci PASSED
   [testng] PASSED: testCountLoci</code class="pre_md"></pre>
<h3>Generating initial MD5s</h3>
<p>If you don't know the MD5s yet you can run the command yourself on the command line and then MD5s the outputs yourself, or you can set the MD5s in your test to <code>""</code> and run the pipeline.</p>
<p>When the MD5s are blank as in:</p>
<pre><code class="pre_md">spec.fileMD5s += testOut -&gt; ""</code class="pre_md"></pre>
<p>You run: </p>
<pre><code class="pre_md">ant pipelinetest -Dsingle=ExampleCountLociPipelineTest -Dpipeline.run=run</code class="pre_md"></pre>
<p>And the output will look like:</p>
<pre><code class="pre_md">   [testng] --------------------------------------------------------------------------------
   [testng] Executing test countloci with Queue arguments: -S scala/qscript/examples/ExampleCountLoci.scala -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/small_bam_for_countloci.bam -o count.out -bsub -l WARN -tempDir pipelinetests/countloci/temp/ -runDir pipelinetests/countloci/run/ -jobQueue hour -run
   [testng] ##### MD5 file is up to date: integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] PARAMETERIZATION[countloci]: file pipelinetests/countloci/run/count.out has md5 = 67823e4722495eb10a5e4c42c267b3a6, stated expectation is , equal? = false
   [testng]   =&gt; countloci PASSED
   [testng] PASSED: testCountLoci</code class="pre_md"></pre>
<h3>Checking MD5s</h3>
<p>When a pipeline test fails due to an MD5 mismatch you can use the MD5 database to diff the results.</p>
<pre><code class="pre_md">   [testng] --------------------------------------------------------------------------------
   [testng] Executing test countloci with Queue arguments: -S scala/qscript/examples/ExampleCountLoci.scala -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/small_bam_for_countloci.bam -o count.out -bsub -l WARN -tempDir pipelinetests/countloci/temp/ -runDir pipelinetests/countloci/run/ -jobQueue hour -run
   [testng] ##### Updating MD5 file: integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] Checking MD5 for pipelinetests/countloci/run/count.out [calculated=67823e4722495eb10a5e4c42c267b3a6, expected=67823e4722495eb10a5e0000deadbeef]
   [testng] ##### Test countloci is going fail #####
   [testng] ##### Path to expected   file (MD5=67823e4722495eb10a5e0000deadbeef): integrationtests/67823e4722495eb10a5e0000deadbeef.integrationtest
   [testng] ##### Path to calculated file (MD5=67823e4722495eb10a5e4c42c267b3a6): integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] ##### Diff command: diff integrationtests/67823e4722495eb10a5e0000deadbeef.integrationtest integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] FAILED: testCountLoci
   [testng] java.lang.AssertionError: 1 of 1 MD5s did not match.</code class="pre_md"></pre>
<p>If you need to examine a number of MD5s which may have changed you can briefly shut off MD5 mismatch failures by setting <code>parameterize = true</code>.</p>
<pre><code class="pre_md">spec.parameterize = true
spec.fileMD5s += testOut -&gt; "67823e4722495eb10a5e4c42c267b3a6"</code class="pre_md"></pre>
<p>For this run:  </p>
<pre><code class="pre_md">ant pipelinetest -Dsingle=ExampleCountLociPipelineTest -Dpipeline.run=run</code class="pre_md"></pre>
<p>If there's a match the output will resemble:</p>
<pre><code class="pre_md">   [testng] --------------------------------------------------------------------------------
   [testng] Executing test countloci with Queue arguments: -S scala/qscript/examples/ExampleCountLoci.scala -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/small_bam_for_countloci.bam -o count.out -bsub -l WARN -tempDir pipelinetests/countloci/temp/ -runDir pipelinetests/countloci/run/ -jobQueue hour -run
   [testng] ##### MD5 file is up to date: integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] PARAMETERIZATION[countloci]: file pipelinetests/countloci/run/count.out has md5 = 67823e4722495eb10a5e4c42c267b3a6, stated expectation is 67823e4722495eb10a5e4c42c267b3a6, equal? = true
   [testng]   =&gt; countloci PASSED
   [testng] PASSED: testCountLoci</code class="pre_md"></pre>
<p>While for a mismatch it will look like this:</p>
<pre><code class="pre_md">   [testng] --------------------------------------------------------------------------------
   [testng] Executing test countloci with Queue arguments: -S scala/qscript/examples/ExampleCountLoci.scala -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/small_bam_for_countloci.bam -o count.out -bsub -l WARN -tempDir pipelinetests/countloci/temp/ -runDir pipelinetests/countloci/run/ -jobQueue hour -run
   [testng] ##### MD5 file is up to date: integrationtests/67823e4722495eb10a5e4c42c267b3a6.integrationtest
   [testng] PARAMETERIZATION[countloci]: file pipelinetests/countloci/run/count.out has md5 = 67823e4722495eb10a5e4c42c267b3a6, stated expectation is 67823e4722495eb10a5e0000deadbeef, equal? = false
   [testng]   =&gt; countloci PASSED
   [testng] PASSED: testCountLoci</code class="pre_md"></pre>
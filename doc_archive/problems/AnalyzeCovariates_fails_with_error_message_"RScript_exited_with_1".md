## AnalyzeCovariates fails with error message "RScript exited with 1"

http://gatkforums.broadinstitute.org/gatk/discussion/4294/analyzecovariates-fails-with-error-message-rscript-exited-with-1

<p>When you run AnalyzeCovariates to analyze your BQSR outputs, you may encounter an error starting with this line:</p>
<pre><code class="pre_md">org.broadinstitute.sting.utils.R.RScriptExecutorException: RScript exited with 1. Run with -l DEBUG for more info.</code class="pre_md"></pre>
<p>The main reason why this error often occurs is simple, and so is the solution. The script depends on some external R libraries, so if you don't have them installed, the script fails. To find out what libraries are necessary and how to install them, you can refer to <a href="http://www.broadinstitute.org/gatk/guide/article?id=2899">this tutorial</a>.</p>
<p>One other common issue is that the version of ggplot2 you have installed is very recent and is not compatible with the BQSR script. If so, download <a href="https://github.com/broadgsa/gatk/blob/6ba57d05eb20101517d8888f999d6d1f564d2aeb/public/gatk-tools-public/src/main/resources/org/broadinstitute/gatk/utils/recalibration/BQSR.R">this Rscript file</a> and use it to generate the plots manually according to the instructions below.</p>
<p>If you have already checked that you have all the necessary libraries installed, you'll need to run the script manually in order to find out what is wrong. To new users, this can seem complicated, but it only takes these 3 simple steps to do it! </p>
<h3>1. Re-run AnalyzeCovariates with these additional parameters:</h3>
<ul>
<li><code>-l DEBUG</code> (that's a lowercase L, not an uppercase i, to be clear) and</li>
<li><code>-csv my-report.csv</code> (where you can call the .csv file anything; this is so the intermediate csv file will be saved).</li>
</ul>
<h3>2. Identify the lines in the log output that says what parameters the RScript is given.</h3>
<p>The snippet below shows you the components of the R script command line that AnalyzeCovariates uses.</p>
<pre><code class="pre_md">INFO  18:04:55,355 AnalyzeCovariates - Generating plots file 'RTest.pdf' 
DEBUG 18:04:55,672 RecalUtils - R command line: Rscript (resource)org/broadinstitute/gatk/utils/recalibration/BQSR.R /Users/schandra/BQSR_Testing/RTest.csv /Users/schandra/BQSR_Testing/RTest.recal /Users/schandra/BQSR_Testing/RTest.pdf 
DEBUG 18:04:55,687 RScriptExecutor - Executing: 
DEBUG 18:04:55,688 RScriptExecutor -   Rscript 
DEBUG 18:04:55,688 RScriptExecutor -   -e 
DEBUG 18:04:55,688 RScriptExecutor -   tempLibDir = '/var/folders/j9/5qgr3mvj0590pd2yb9hwc15454pxz0/T/Rlib.2085451458391709180';source('/var/folders/j9/5qgr3mvj0590pd2yb9hwc15454pxz0/T/BQSR.761775214345441497.R'); 
DEBUG 18:04:55,689 RScriptExecutor -   /Users/schandra/BQSR_Testing/RTest.csv 
DEBUG 18:04:55,689 RScriptExecutor -   /Users/schandra/BQSR_Testing/RTest.recal 
DEBUG 18:04:55,689 RScriptExecutor -   /Users/schandra/BQSR_Testing/RTest.pdf </code class="pre_md"></pre>
<p>So, your full command line will be:</p>
<pre><code class="pre_md">RScript BQSR.R RTest.csv RTest.recal RTest.pdf</code class="pre_md"></pre>
<p><strong>Please note:</strong></p>
<ul>
<li>BQSR.R is the name of the script you want to run. It can be found [here]( <a href="https://github.com/broadgsa/gatk/blob/6ba57d05eb20101517d8888f999d6d1f564d2aeb/public/gatk-tools-public/src/main/resources/org/broadinstitute/gatk/utils/recalibration/BQSR.R">https://github.com/broadgsa/gatk/blob/6ba57d05eb20101517d8888f999d6d1f564d2aeb/public/gatk-tools-public/src/main/resources/org/broadinstitute/gatk/utils/recalibration/BQSR.R</a>)  </li>
<li>RTest.csv is the name of the original csv file output from AnalyzeCovariates.  </li>
<li>RTest.recal is your original recalibration file.   </li>
<li>RTest.pdf is the output pdf file; you can name it whatever you want.  </li>
</ul>
<h3>3. Run the script manually with the above arguments.</h3>
<p>For new users, the easiest way to do this is to do it from within an IDE program like RStudio. Or, you can start up R at the command line and run it that way, whatever you are comfortable with. </p>
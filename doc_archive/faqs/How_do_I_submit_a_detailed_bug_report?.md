## How do I submit a detailed bug report?

http://gatkforums.broadinstitute.org/gatk/discussion/1894/how-do-i-submit-a-detailed-bug-report

<p><strong><em>Note: only do this if you have been explicitly asked to do so.</em></strong></p>
<h3>Scenario:</h3>
<p>You posted a question about a problem you had with GATK tools, we answered that we think it's a bug, and we asked you to submit a detailed bug report. </p>
<h3>Here's what you need to provide:</h3>
<ul>
<li>The exact command line that you used when you had the problem (in a text file)</li>
<li>The full log output (program output in the console) from the start of the run to the end or error message (in a text file)</li>
<li>A snippet of the BAM file if applicable and the index (.bai) file associated with it</li>
<li>If a non-standard reference (i.e. not available in our resource bundle) was used, we need the .fasta, .fai, and .dict files for the reference</li>
<li>Any other relevant files such as recalibration plots</li>
</ul>
<p>A snippet file is a slice of the original BAM file which contains the problematic region and is sufficient to reproduce the error. We need it in order to reproduce the problem on our end, which is the first necessary step to finding and fixing the bug. We ask you to provide this as a snippet rather than the full file so that you don't have to upload (and we don't have to process) huge giga-scale files.  </p>
<h3>Here's how you create a snippet file:</h3>
<ul>
<li>Look at the error message and see if it cites a specific position where the error occurred</li>
<li>If not, identify what region caused the problem by running with <code>-L</code> argument and progressively narrowing down the interval</li>
<li>Once you have the region, use PrintReads with <code>-L</code> to write the problematic region (with 500 bp padding on either side) to a new file -- this is your snippet file.</li>
<li>Test your command line on this snippet  file to make sure you can still reproduce the error on it. </li>
</ul>
<h3>And finally, here's how you send us the files:</h3>
<ul>
<li>Put all those files into a <code>.zip</code> or <code>.tar.gz</code> archive </li>
<li>
<p>Upload them onto our FTP server with the following credentials:</p>
<pre><code>location: ftp.broadinstitute.org
username: gsapubftp
password: 5WvQWSfi</code></pre>
</li>
<li>Post in the original discussion thread that you have done this</li>
<li>Be sure to tell us the name of your archive file!</li>
</ul>
<h3>We will get back to you --hopefully with a bug fix!-- as soon as we can.</h3>
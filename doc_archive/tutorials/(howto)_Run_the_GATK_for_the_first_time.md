## (howto) Run the GATK for the first time

http://gatkforums.broadinstitute.org/gatk/discussion/1209/howto-run-the-gatk-for-the-first-time

<h4>NOTICE:</h4>
<p>This tutorial is slightly out of date so the output is a little different. We'll update this soon, but in the meantime, don't freak out if you get a result that reads something like </p>
<pre><code class="pre_md">INFO 18:32:38,826 CountReads - CountReads counted 33 reads in the traversal </code class="pre_md"></pre>
<p>instead of </p>
<pre><code class="pre_md">INFO  16:17:46,061 Walker - [REDUCE RESULT] Traversal result is: 33 </code class="pre_md"></pre>
<p>You're doing the right thing and getting the right result.</p>
<p>And of course, in doubt, just post a comment on this article; we're here to answer your questions. </p>
<hr />
<h4>Objective</h4>
<p>Run a basic analysis command on example data.</p>
<h4>Prerequisites</h4>
<ul>
<li>Successfully completed <a href="http://gatkforums.broadinstitute.org/discussion/1200/how-to-test-your-gatk-installation">&quot;How to test your GATK installation&quot;</a></li>
<li>Familiarity with <a href="http://gatkforums.broadinstitute.org/discussion/1204/what-input-files-does-the-gatk-accept">&quot;Input files for the GATK&quot;</a></li>
<li><a href="http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it">GATK resource bundle</a> downloaded </li>
</ul>
<h4>Steps</h4>
<ol>
<li>Invoke the GATK CountReads command</li>
<li>Further exercises</li>
</ol>
<hr />
<h3>1. Invoke the GATK CountReads command</h3>
<p>A very simple analysis that you can do with the GATK is getting a count of the reads in a BAM file. The GATK is capable of much more powerful analyses, but this is a good starting example because there are very few things that can go wrong.</p>
<p>So we are going to count the reads in the file <code>exampleBAM.bam</code>, which you can find in the <a href="http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it">GATK resource bundle</a> along with its associated index (same file name with <code>.bai</code> extension), as well as the example reference <code>exampleFASTA.fasta</code> and its associated index (same file name with <code>.fai</code> extension) and dictionary (same file name with <code>.dict</code> extension). Copy them to your working directory so that your directory contents look like this:</p>
<pre><code class="pre_md">[bm4dd-56b:~/codespace/gatk/sandbox] vdauwera% ls -la
drwxr-xr-x  9 vdauwera  CHARLES\Domain Users     306 Jul 25 16:29 .
drwxr-xr-x@ 6 vdauwera  CHARLES\Domain Users     204 Jul 25 15:31 ..
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users    3635 Apr 10 07:39 exampleBAM.bam
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users     232 Apr 10 07:39 exampleBAM.bam.bai
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users     148 Apr 10 07:39 exampleFASTA.dict
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users  101673 Apr 10 07:39 exampleFASTA.fasta
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users      20 Apr 10 07:39 exampleFASTA.fasta.fai</code class="pre_md"></pre>
<h4>Action</h4>
<p>Type the following command:</p>
<pre><code class="pre_md">java -jar &lt;path to GenomeAnalysisTK.jar&gt; -T CountReads -R exampleFASTA.fasta -I exampleBAM.bam </code class="pre_md"></pre>
<p>where <code>-T CountReads</code> specifies which analysis tool we want to use, <code>-R exampleFASTA.fasta</code> specifies the reference sequence, and <code>-I exampleBAM.bam</code> specifies the file of aligned reads we want to analyze.</p>
<p>For any analysis that you want to run on a set of aligned reads, you will <strong>always</strong> need to use at least these three arguments: </p>
<ul>
<li><code>-T</code> for the <strong>t</strong>ool name, which specifices the corresponding analysis</li>
<li><code>-R</code> for the <strong>r</strong>eference sequence file</li>
<li><code>-I</code> for the <strong>i</strong>nput BAM file of aligned reads</li>
</ul>
<p>They don't have to be in that order in your command, but this way you can remember that you need them if you <strong>TRI</strong>...</p>
<h4>Expected Result</h4>
<p>After a few seconds you should see output that looks like to this:</p>
<pre><code class="pre_md">INFO  16:17:45,945 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:17:45,946 HelpFormatter - The Genome Analysis Toolkit (GATK) v2.0-22-g40f97eb, Compiled 2012/07/25 15:29:41 
INFO  16:17:45,947 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  16:17:45,947 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  16:17:45,947 HelpFormatter - Program Args: -T CountReads -R exampleFASTA.fasta -I exampleBAM.bam 
INFO  16:17:45,947 HelpFormatter - Date/Time: 2012/07/25 16:17:45 
INFO  16:17:45,947 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:17:45,948 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:17:45,950 GenomeAnalysisEngine - Strictness is SILENT 
INFO  16:17:45,982 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  16:17:45,993 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.01 
INFO  16:17:46,060 TraversalEngine - [INITIALIZATION COMPLETE; TRAVERSAL STARTING] 
INFO  16:17:46,060 TraversalEngine -        Location processed.reads  runtime per.1M.reads completed total.runtime remaining 
INFO  16:17:46,061 Walker - [REDUCE RESULT] Traversal result is: 33 
INFO  16:17:46,061 TraversalEngine - Total runtime 0.00 secs, 0.00 min, 0.00 hours 
INFO  16:17:46,100 TraversalEngine - 0 reads were filtered out during traversal out of 33 total (0.00%) 
INFO  16:17:46,729 GATKRunReport - Uploaded run statistics report to AWS S3 </code class="pre_md"></pre>
<p>Depending on the GATK release, you may see slightly different information output, but you know everything is running correctly if you see the line:</p>
<pre><code class="pre_md">INFO  21:53:04,556 Walker - [REDUCE RESULT] Traversal result is: 33 </code class="pre_md"></pre>
<p>somewhere in your output.</p>
<p>If you don't see this, check your spelling (GATK commands are case-sensitive), check that the files are in your working directory, and if necessary, re-check that the GATK is properly installed.</p>
<p>If you do see this output, congratulations! You just successfully ran you first GATK analysis! </p>
<p>Basically the output you see means that the CountReadsWalker (which you invoked with the command line option <code>-T CountReads</code>) counted 33 reads in the <code>exampleBAM.bam</code> file, which is exactly what we expect to see. </p>
<p><strong>Wait, what is this <em>walker</em> thing?</strong></p>
<p>In the GATK jargon, we call the tools <em>walkers</em> because the way they work is that they <em>walk</em> through the dataset --either along the reference sequence (LocusWalkers), or down the list of reads in the BAM file (ReadWalkers)-- collecting the requested information along the way. </p>
<hr />
<h3>2. Further Exercises</h3>
<p>Now that you're rocking the read counts, you can start to expand your use of the GATK command line.</p>
<p>Let's say you don't care about counting reads anymore; now you want to know the number of loci (positions on the genome) that are covered by one or more reads. The name of the tool, or walker, that does this is <code>CountLoci</code>. Since the structure of the GATK command is basically always the same, you can simply switch the tool name, right?</p>
<h4>Action</h4>
<p>Instead of this command, which we used earlier:</p>
<pre><code class="pre_md">java -jar &lt;path to GenomeAnalysisTK.jar&gt; -T CountReads -R exampleFASTA.fasta -I exampleBAM.bam </code class="pre_md"></pre>
<p>this time you type this:</p>
<pre><code class="pre_md">java -jar &lt;path to GenomeAnalysisTK.jar&gt; -T CountLoci -R exampleFASTA.fasta -I exampleBAM.bam </code class="pre_md"></pre>
<p>See the difference?</p>
<h4>Result</h4>
<p>You should see something like this output:</p>
<pre><code class="pre_md">INFO  16:18:26,183 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:18:26,185 HelpFormatter - The Genome Analysis Toolkit (GATK) v2.0-22-g40f97eb, Compiled 2012/07/25 15:29:41 
INFO  16:18:26,185 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  16:18:26,185 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  16:18:26,186 HelpFormatter - Program Args: -T CountLoci -R exampleFASTA.fasta -I exampleBAM.bam 
INFO  16:18:26,186 HelpFormatter - Date/Time: 2012/07/25 16:18:26 
INFO  16:18:26,186 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:18:26,186 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:18:26,189 GenomeAnalysisEngine - Strictness is SILENT 
INFO  16:18:26,222 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  16:18:26,233 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.01 
INFO  16:18:26,351 TraversalEngine - [INITIALIZATION COMPLETE; TRAVERSAL STARTING] 
INFO  16:18:26,351 TraversalEngine -        Location processed.sites  runtime per.1M.sites completed total.runtime remaining 
2052
INFO  16:18:26,411 TraversalEngine - Total runtime 0.08 secs, 0.00 min, 0.00 hours 
INFO  16:18:26,450 TraversalEngine - 0 reads were filtered out during traversal out of 33 total (0.00%) 
INFO  16:18:27,124 GATKRunReport - Uploaded run statistics report to AWS S3 </code class="pre_md"></pre>
<p>Great! But wait -- where's the result? Last time the result was given on this line:</p>
<pre><code class="pre_md">INFO  21:53:04,556 Walker - [REDUCE RESULT] Traversal result is: 33 </code class="pre_md"></pre>
<p><strong>But this time there is no line that says <code>[REDUCE RESULT]</code>! Is something wrong?</strong></p>
<p>Not really. The program ran just fine -- but we forgot to give it an output file name. You see, the <code>CountLoci</code> walker is set up to output the result of its calculations to a text file, unlike <code>CountReads</code>, which is perfectly happy to output its result to the terminal screen. </p>
<h4>Action</h4>
<p>So we repeat the command, but this time we specify an output file, like this:</p>
<pre><code class="pre_md">java -jar &lt;path to GenomeAnalysisTK.jar&gt; -T CountLoci -R exampleFASTA.fasta -I exampleBAM.bam -o output.txt</code class="pre_md"></pre>
<p>where <code>-o</code> (lowercase o, not zero) is used to specify the output.  </p>
<h4>Result</h4>
<p>You should get essentially the same output on the terminal screen as previously (but notice the difference in the line that contains <code>Program Args</code> -- the new argument is included):</p>
<pre><code class="pre_md">INFO  16:29:15,451 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:29:15,453 HelpFormatter - The Genome Analysis Toolkit (GATK) v2.0-22-g40f97eb, Compiled 2012/07/25 15:29:41 
INFO  16:29:15,453 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  16:29:15,453 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  16:29:15,453 HelpFormatter - Program Args: -T CountLoci -R exampleFASTA.fasta -I exampleBAM.bam -o output.txt 
INFO  16:29:15,454 HelpFormatter - Date/Time: 2012/07/25 16:29:15 
INFO  16:29:15,454 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:29:15,454 HelpFormatter - --------------------------------------------------------------------------------- 
INFO  16:29:15,457 GenomeAnalysisEngine - Strictness is SILENT 
INFO  16:29:15,488 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  16:29:15,499 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.01 
INFO  16:29:15,618 TraversalEngine - [INITIALIZATION COMPLETE; TRAVERSAL STARTING] 
INFO  16:29:15,618 TraversalEngine -        Location processed.sites  runtime per.1M.sites completed total.runtime remaining 
INFO  16:29:15,679 TraversalEngine - Total runtime 0.08 secs, 0.00 min, 0.00 hours 
INFO  16:29:15,718 TraversalEngine - 0 reads were filtered out during traversal out of 33 total (0.00%) 
INFO  16:29:16,712 GATKRunReport - Uploaded run statistics report to AWS S3 </code class="pre_md"></pre>
<p>This time however, if we look inside the working directory, there is a newly created file there called <code>output.txt</code>. </p>
<pre><code class="pre_md">[bm4dd-56b:~/codespace/gatk/sandbox] vdauwera% ls -la
drwxr-xr-x  9 vdauwera  CHARLES\Domain Users     306 Jul 25 16:29 .
drwxr-xr-x@ 6 vdauwera  CHARLES\Domain Users     204 Jul 25 15:31 ..
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users    3635 Apr 10 07:39 exampleBAM.bam
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users     232 Apr 10 07:39 exampleBAM.bam.bai
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users     148 Apr 10 07:39 exampleFASTA.dict
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users  101673 Apr 10 07:39 exampleFASTA.fasta
-rw-r--r--@ 1 vdauwera  CHARLES\Domain Users      20 Apr 10 07:39 exampleFASTA.fasta.fai
-rw-r--r--  1 vdauwera  CHARLES\Domain Users       5 Jul 25 16:29 output.txt</code class="pre_md"></pre>
<p>This file contains the result of the analysis:</p>
<pre><code class="pre_md">[bm4dd-56b:~/codespace/gatk/sandbox] vdauwera% cat output.txt 
2052</code class="pre_md"></pre>
<p>This means that there are 2052 loci in the reference sequence that are covered by at least one or more reads in the BAM file.  </p>
<h4>Discussion</h4>
<p>Okay then, but why not show the full, correct command in the first place? Because this was a good opportunity for you to learn a few of the caveats of the GATK command system, which may save you a lot of frustration later on.</p>
<p>Beyond the common basic arguments that almost all GATK walkers require, most of them also have specific requirements or options that are important to how they work. You should always check what are the specific arguments that are required, recommended and/or optional for the walker you want to use before starting an analysis. </p>
<p>Fortunately the GATK is set up to complain (<em>i.e.</em> terminate with an error message) if you try to run it without specifying a required argument. For example, if you try to run this:</p>
<pre><code class="pre_md">java -jar &lt;path to GenomeAnalysisTK.jar&gt; -T CountLoci -R exampleFASTA.fasta</code class="pre_md"></pre>
<p>the GATK will spit out a wall of text, including the basic usage guide that you can invoke with the <code>--help</code> option, and more importantly, the following error message: </p>
<pre><code class="pre_md">##### ERROR ------------------------------------------------------------------------------------------
##### ERROR A USER ERROR has occurred (version 2.0-22-g40f97eb): 
##### ERROR The invalid arguments or inputs must be corrected before the GATK can proceed
##### ERROR Please do not post this error to the GATK forum
##### ERROR
##### ERROR See the documentation (rerun with -h) for this tool to view allowable command-line arguments.
##### ERROR Visit our website and forum for extensive documentation and answers to 
##### ERROR commonly asked questions http://www.broadinstitute.org/gatk
##### ERROR
##### ERROR MESSAGE: Walker requires reads but none were provided.
##### ERROR ------------------------------------------------------------------------------------------</code class="pre_md"></pre>
<p>You see the line that says <code>ERROR MESSAGE: Walker requires reads but none were provided</code>? This tells you exactly what was wrong with your command. </p>
<p>So the GATK will not run if a walker does not have all the required inputs. That's a good thing! But in the case of our first attempt at running <code>CountLoci</code>, the <code>-o</code> argument is not required by the GATK to run -- it's just highly desirable if you actually want the result of the analysis! </p>
<p>There will be many other cases of walkers with arguments that are not strictly required, but highly desirable if you want the results to be meaningful. </p>
<p>So, at the risk of getting repetitive, <strong>always read the <a href="http://www.broadinstitute.org/gatk/gatkdocs/">documentation</a> of each walker that you want to use!</strong> </p>
## What is the GATKReport file format?

http://gatkforums.broadinstitute.org/gatk/discussion/1244/what-is-the-gatkreport-file-format

<p>A GATKReport is simply a text document that contains well-formatted, easy to read representation of some tabular data. Many GATK tools output their results as GATKReports, so it's important to understand how they are formatted and how you can use them in further analyses.</p>
<p>Here's a simple example:</p>
<pre><code class="pre_md">#:GATKReport.v1.0:2
#:GATKTable:true:2:9:%.18E:%.15f:;
#:GATKTable:ErrorRatePerCycle:The error rate per sequenced position in the reads
cycle  errorrate.61PA8.7         qualavg.61PA8.7                                         
0      7.451835696110506E-3      25.474613284804366                                      
1      2.362777171937477E-3      29.844949954504095                                      
2      9.087604507451836E-4      32.875909752547310
3      5.452562704471102E-4      34.498999090081895                                      
4      9.087604507451836E-4      35.148316651501370                                       
5      5.452562704471102E-4      36.072234352256190                                       
6      5.452562704471102E-4      36.121724890829700                                        
7      5.452562704471102E-4      36.191048034934500                                        
8      5.452562704471102E-4      36.003457059679770                                       

#:GATKTable:false:2:3:%s:%c:;
#:GATKTable:TableName:Description
key    column
1:1000  T 
1:1001  A 
1:1002  C </code class="pre_md"></pre>
<p>This report contains two individual GATK report tables. Every table begins with a header for its metadata and then a header for its name and description. The next row contains the column names followed by the data. </p>
<p>We provide an R library called <code>gsalib</code> that allows you to load GATKReport files into R for further analysis. Here are four simple steps to getting <code>gsalib</code>, installing it and loading a report. </p>
<h4>1. Start R (or open RStudio)</h4>
<pre><code class="pre_md">$ R

R version 2.11.0 (2010-04-22)
Copyright (C) 2010 The R Foundation for Statistical Computing
ISBN 3-900051-07-0

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.</code class="pre_md"></pre>
<h4>2. Get the <code>gsalib</code> library from CRAN</h4>
<p>The <code>gsalib</code> library is available on the <a href="http://cran.r-project.org/">Comprehensive R Archive Network</a>, so you can just do:</p>
<pre><code class="pre_md">&gt; install.packages("gsalib") </code class="pre_md"></pre>
<p>From within R (we use RStudio for convenience).</p>
<p>In some cases you need to explicitly tell R where to find the library; you can do this as follows:</p>
<pre><code class="pre_md">$ cat .Rprofile 
.libPaths("/path/to/Sting/R/")</code class="pre_md"></pre>
<h4>3. Load the gsalib library</h4>
<pre><code class="pre_md">&gt; library(gsalib)</code class="pre_md"></pre>
<h4>4. Finally, load the GATKReport file and have fun</h4>
<pre><code class="pre_md">&gt; d = gsa.read.gatkreport("/path/to/my.gatkreport")
&gt; summary(d)
              Length Class      Mode
CountVariants 27     data.frame list
CompOverlap   13     data.frame list</code class="pre_md"></pre>
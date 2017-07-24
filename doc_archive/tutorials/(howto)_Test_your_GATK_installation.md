## (howto) Test your GATK installation

http://gatkforums.broadinstitute.org/gatk/discussion/1200/howto-test-your-gatk-installation

<h4>Objective</h4>
<p>Test that the GATK is correctly installed, and that the supporting tools like Java are in your path.</p>
<h4>Prerequisites</h4>
<ul>
<li>Basic familiarity with the command-line environment</li>
<li>Understand what is a PATH variable</li>
<li>GATK downloaded and placed on path</li>
</ul>
<h4>Steps</h4>
<ol>
<li>Invoke the GATK usage/help message</li>
<li>Troubleshooting</li>
</ol>
<hr />
<h3>1. Invoke the GATK usage/help message</h3>
<p>The command we're going to run is a very simple command that asks the GATK to print out a list of available command-line arguments and options. It is so simple that it will ALWAYS work if your GATK package is installed correctly.</p>
<p>Note that this command is also helpful when you're trying to remember something like the right spelling or short name for an argument and for whatever reason you don't have access to the web-based documentation.  </p>
<h4>Action</h4>
<p>Type the following command:</p>
<pre><code class="pre_md">java -jar &lt;path to GenomeAnalysisTK.jar&gt; --help</code class="pre_md"></pre>
<p>replacing the <code>&lt;path to GenomeAnalysisTK.jar&gt;</code> bit with the path you have set up in your command-line environment.</p>
<h4>Expected Result</h4>
<p>You should see usage output similar to the following:</p>
<pre><code class="pre_md">usage: java -jar GenomeAnalysisTK.jar -T &lt;analysis_type&gt; [-I &lt;input_file&gt;] [-L 
        &lt;intervals&gt;] [-R &lt;reference_sequence&gt;] [-B &lt;rodBind&gt;] [-D &lt;DBSNP&gt;] [-H 
        &lt;hapmap&gt;] [-hc &lt;hapmap_chip&gt;] [-o &lt;out&gt;] [-e &lt;err&gt;] [-oe &lt;outerr&gt;] [-A] [-M 
        &lt;maximum_reads&gt;] [-sort &lt;sort_on_the_fly&gt;] [-compress &lt;bam_compression&gt;] [-fmq0] [-dfrac 
        &lt;downsample_to_fraction&gt;] [-dcov &lt;downsample_to_coverage&gt;] [-S 
        &lt;validation_strictness&gt;] [-U] [-P] [-dt] [-tblw] [-nt &lt;numthreads&gt;] [-l 
        &lt;logging_level&gt;] [-log &lt;log_to_file&gt;] [-quiet] [-debug] [-h]
-T,--analysis_type &lt;analysis_type&gt;                     Type of analysis to run
-I,--input_file &lt;input_file&gt;                           SAM or BAM file(s)
-L,--intervals &lt;intervals&gt;                             A list of genomic intervals over which 
                                                       to operate. Can be explicitly specified 
                                                       on the command line or in a file.
-R,--reference_sequence &lt;reference_sequence&gt;           Reference sequence file
-B,--rodBind &lt;rodBind&gt;                                 Bindings for reference-ordered data, in 
                                                       the form &lt;name&gt;,&lt;type&gt;,&lt;file&gt;
-D,--DBSNP &lt;DBSNP&gt;                                     DBSNP file
-H,--hapmap &lt;hapmap&gt;                                   Hapmap file
-hc,--hapmap_chip &lt;hapmap_chip&gt;                        Hapmap chip file
-o,--out &lt;out&gt;                                         An output file presented to the walker. 
                                                       Will overwrite contents if file exists.
-e,--err &lt;err&gt;                                         An error output file presented to the 
                                                       walker. Will overwrite contents if file 
                                                       exists.
-oe,--outerr &lt;outerr&gt;                                  A joint file for 'normal' and error 
                                                       output presented to the walker. Will 
                                                       overwrite contents if file exists.

...</code class="pre_md"></pre>
<p>If you see this message, your GATK installation is ok. You're good to go! If you don't see this message, and instead get an error message, proceed to the next section on troubleshooting.  </p>
<hr />
<h3>2. Troubleshooting</h3>
<p>Let's try to figure out what's not working.  </p>
<h4>Action</h4>
<p>First, make sure that your Java version is at least 1.7, by typing the following command:</p>
<pre><code class="pre_md">java -version</code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>You should see something similar to the following text:</p>
<pre><code class="pre_md">java version "1.7.0_12"
Java(TM) SE Runtime Environment (build 1.7.0_12-b04)
Java HotSpot(TM) 64-Bit Server VM (build 11.2-b01, mixed mode)  </code class="pre_md"></pre>
<h4>Remedial actions</h4>
<p>If the version is less then 1.7, install the newest version of Java onto the system. If you instead see something like </p>
<pre><code class="pre_md">java: Command not found  </code class="pre_md"></pre>
<p>make sure that java is installed on your machine, and that your PATH variable contains the path to the java executables. </p>
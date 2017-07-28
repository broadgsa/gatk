## (howto) Install software for GATK workshops

http://gatkforums.broadinstitute.org/gatk/discussion/7098/howto-install-software-for-gatk-workshops

<h3>Objective</h3>
<p>Install all software packages required to attend a GATK workshop. </p>
<h3>Prerequisites</h3>
<p>To follow these instructions, you will need to have a basic understanding of the meaning of the following words and command-line operations. If you are unfamiliar with any of the following, you should consult a more experienced colleague or your system administrator if you have one. There are also many good online tutorials you can use to learn the necessary notions.</p>
<ul>
<li>Basic Unix environment commands </li>
<li>Binary / Executable </li>
<li>Adding a binary to your path (optional)</li>
<li>Command-line shell, terminal or console </li>
<li>Software library</li>
</ul>
<h3>Platform requirements</h3>
<p>GATK is supported on all flavors of reasonably recent Linux/Unix and MacOS X systems, but <strong>NOT on Windows</strong>. The analyses we run in workshops are designed to run quickly and on small datasets, so should not require more than 2G of RAM. For file storage, plan on 10G of space (but I would be shocked if we get to half of that).</p>
<p>The current version of GATK requires Java Runtime Environment version 1.8. All Linux/Unix and MacOS X systems should have a JRE pre-installed, but the version may vary. To test your Java version, run the following command in the shell: </p>
<pre><code class="pre_md">java -version </code class="pre_md"></pre>
<p>This should return a message along the lines of ”java version 1.8.0_65” as well as some details on the Runtime Environment (JRE) and Virtual Machine (VM). If you have a version other than 1.8.x, be aware that you may run into trouble with some of the more advanced features of the Picard and GATK tools. The simplest solution is to install an additional JRE and specify which you want to use at the command-line. To find out how to do so, you should seek help from your system administrator and read <a href="https://www.broadinstitute.org/gatk/guide/article?id=6841">this article</a>. </p>
<h3>Software packages</h3>
<ol>
<li>Picard</li>
<li>Genome Analysis Toolkit (GATK) </li>
<li>IGV  </li>
<li>RStudio IDE and R libraries ggplot2 and gsalib  </li>
<li>Samtools</li>
<li>RTG Tools</li>
</ol>
<hr />
<h3>1. Picard</h3>
<p>Read the overview of the Picard software on the <a href="http://broadinstitute.github.io/picard/">Picard project homepage</a>, then download the <a href="https://github.com/broadinstitute/picard/releases/">latest version</a> (currently 2.4.1) of the package containing the pre-compiled program file (the picard-tools-2.x.y.zip file). </p>
<ul>
<li>Installation</li>
</ul>
<p>Unpack the zip file using:   </p>
<pre><code class="pre_md">tar xjf picard-tools-2.4.1.zip </code class="pre_md"></pre>
<p>This will produce a directory called <code>picard-tools-2.4.1</code> containing the Picard jar files. Picard tools are distributed as a pre-compiled Java executable (jar file) so there is no need to compile them. </p>
<p>Note that it is not possible to add jar files to your path to make the tools available on the command line; you have to specify the full path to the jar file in your java command, which would look like this: </p>
<pre><code class="pre_md">java -jar ~/my_tools/jars/picard.jar &lt;Toolname&gt; [options]</code class="pre_md"></pre>
<p><em>This syntax will be explained in a little more detail further below.</em></p>
<p>However, you can set up a shortcut called an &quot;environment variable&quot; in your shell profile configuration to make this easier. The idea is that you create a variable that tells your system where to find a given jar, like this:</p>
<pre><code class="pre_md">PICARD = "~/my_tools/jars/picard.jar"</code class="pre_md"></pre>
<p>So then when you want to run a Picard tool, you just need to call the jar by its shortcut, like this:</p>
<pre><code class="pre_md">java -jar $PICARD &lt;Toolname&gt; [options]</code class="pre_md"></pre>
<p>The exact way to set this up depends on what shell you're using and how your environment is configured. We like <a href="https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps">this overview and tutorial</a> which explains how it all works; but if you are new to the command line environment and you find this too much too deal with, we recommend asking for help from your institution's IT support group.</p>
<p>This completes the installation process. </p>
<ul>
<li>Testing</li>
</ul>
<p>Open a shell and run:  </p>
<pre><code class="pre_md">java -jar picard.jar -h </code class="pre_md"></pre>
<p>This should print out some version and usage information about the <code>AddOrReplaceReadGroups.jar</code> tool. At this point you will have noticed an important difference between BWA and Picard tools. To use BWA, we called on the BWA program and specified which of its internal tools we wanted to apply. To use Picard, we called on Java itself as the main program, then specified which jar file to use, knowing that one jar file = one tool. This applies to all Picard tools; to use them you will always build your command lines like this:   </p>
<pre><code class="pre_md">java -jar picard.jar &lt;ToolName&gt; [options] </code class="pre_md"></pre>
<p>This means you first make the call to Java itself as the main program, then specify the <code>picard.jar</code> file, then specify which tool you want, and finally you pass whatever other arguments (input files, parameters etc.) are needed for the analysis. </p>
<p>Note that the command-line syntax of Picard tools has recently changed from <code>java -jar &lt;ToolName&gt;.jar</code> to <code>java -jar picard.jar &lt;ToolName&gt;</code>. We are using the newer syntax in this document, but some of our other documents may not have been updated yet. If you encounter any documents using the old syntax, let us know and we'll update them accordingly. If you are already using an older version of Picard, either adapt the commands or better, upgrade your version!</p>
<p>Next we will see that GATK tools are called in essentially the same way, although the way the options are specified is a little different. The reasons for how tools in a given software package are organized and invoked are largely due to the preferences of the software developers. They generally do not reflect strict technical requirements, although they can have an effect on speed and efficiency.</p>
<hr />
<h3>2. Genome Analysis Toolkit (GATK)</h3>
<p>Hopefully if you're reading this, you're already acquainted with the <a href="http://www.broadinstitute.org/gatk/about">purpose of the GATK</a>, so go ahead and download the <a href="http://www.broadinstitute.org/gatk/download">latest version of the software package</a>. </p>
<p>In order to access the downloads, you need to register for a free account on the <a href="http://gatkforums.broadinstitute.org/">GATK support forum</a>. You will also need to read and accept the license agreement before downloading the GATK software package. Note that if you intend to use the GATK for commercial purposes, you will need to purchase a license. See the <a href="https://www.broadinstitute.org/gatk/about/#licensing">licensing page</a> for an overview of the commercial licensing conditions. </p>
<ul>
<li>Installation</li>
</ul>
<p>Unpack the tar file using:  </p>
<pre><code class="pre_md">tar xjf GenomeAnalysisTK-3.6-0.tar.bz2 </code class="pre_md"></pre>
<p>This will produce a directory called <code>GenomeAnalysisTK-3.6-0</code> containing the GATK jar file, which is called <code>GenomeAnalysisTK.jar</code>, as well as a directory of example files called <code>resources</code>. GATK tools are distributed as a single pre-compiled Java executable so there is no need to compile them. Just like we discussed for Picard, it's not possible to add the GATK to your path, but you can set up a shortcut to the jar file using environment variables as described above. </p>
<p>This completes the installation process. </p>
<ul>
<li>Testing</li>
</ul>
<p>Open a shell and run: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -h </code class="pre_md"></pre>
<p>This should print out some version and usage information, as well as a list of the tools included in the GATK. As the <strong>Usage</strong> line states, to use GATK you will always build your command lines like this: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T &lt;ToolName&gt; [arguments] </code class="pre_md"></pre>
<p>This means that just like for Picard, you first make the call to Java itself as the main program, then specify the <code>GenomeAnalysisTK.jar</code> file, then specify which tool you want, and finally you pass whatever other arguments (input files, parameters etc.) are needed for the analysis. </p>
<hr />
<h3>3. IGV</h3>
<p>The Integrated Genomics Viewer is a genome browser that allows you to view BAM, VCF and other genomic file information in context. It has a graphical user interface that is very easy to use, and can be downloaded for free (though registration is required) from <a href="http://software.broadinstitute.org/software/igv/">this website</a>. We encourage you to read through IGV's very helpful <a href="http://software.broadinstitute.org/software/igv/UserGuide">user guide</a>, which includes many detailed tutorials that will help you use the program most effectively. </p>
<hr />
<h3>4. RStudio IDE and R libraries ggplot2 and gsalib</h3>
<p>Download the <a href="http://www.rstudio.com/">latest version of RStudio IDE</a>. The webpage should automatically detect what platform you are running on and recommend the version most suitable for your system. </p>
<ul>
<li>Installation</li>
</ul>
<p>Follow the installation instructions provided. Binaries are provided for all major platforms; typically they just need to be placed in your Applications (or Programs) directory. Open RStudio and type the following command in the console window:  </p>
<pre><code class="pre_md">install.packages("ggplot2") </code class="pre_md"></pre>
<p>This will download and install the ggplot2 library as well as any other library packages that ggplot2 depends on for its operation. Note that some users have reported having to install two additional package themselves, called <code>reshape</code> and <code>gplots</code>, which you can do as follows:</p>
<pre><code class="pre_md">install.packages("reshape")
install.packages("gplots")</code class="pre_md"></pre>
<p>Finally, do the same thing to install the gsalib library: </p>
<pre><code class="pre_md">install.packages("gsalib")</code class="pre_md"></pre>
<p>This will download and install the gsalib library.</p>
<hr />
<h3>5. SAMtools</h3>
<p>Read the overview of the SAMtools software on the <a href="http://samtools.sourceforge.net/">SAMtools project homepage</a>, then download the <a href="http://sourceforge.net/projects/samtools/files/">latest version of the software package</a>. </p>
<ul>
<li>Installation</li>
</ul>
<p>Unpack the tar file using:   </p>
<pre><code class="pre_md">tar xvjf samtools-0.1.2.tar.bz2 </code class="pre_md"></pre>
<p>This will produce a directory called <code>samtools-0.1.2</code> containing the files necessary to compile the SAMtools binary. Move to this directory and compile using: </p>
<pre><code class="pre_md">cd samtools-0.1.2 
make </code class="pre_md"></pre>
<p>The compiled binary is called <code>samtools</code>. You should find it within the same folder (<code>samtools-0.1.2</code> in this example). Finally, add the SAMtools binary to your path to make it available on the command line. This completes the installation process. </p>
<ul>
<li>Testing</li>
</ul>
<p>Open a shell and run:  </p>
<pre><code class="pre_md">samtools </code class="pre_md"></pre>
<p>This should print out some version information as well as a list of commands. As the <strong>Usage</strong> line states, to use SAMtools you will always build your command lines like this: </p>
<pre><code class="pre_md">samtools &lt;command&gt; [options] </code class="pre_md"></pre>
<p>This means you first make the call to the binary (<code>samtools</code>), then you specify which command (method) you wish to use (e.g. <code>index</code>) then any options (<em>i.e.</em> arguments such as input files or parameters) used by the program to perform that command. This is a similar convention as used by BWA.</p>
<hr />
<h3>6. RTG Tools</h3>
<p>RTG Tools is a free open-source software package produced by a commercial company called <a href="http://realtimegenomics.com/products/rtg-tools/">Real Time Genomics</a>. This toolkit includes some variant evaluation and plotting tools that we find useful for teaching because they're fairly user-friendly and produce neat interactive plots.</p>
<p>You can download the toolkit from the <a href="http://realtimegenomics.com/products/rtg-tools">RTG website</a>, which provides packages for Linux, MacOS X and Windows. </p>
<ul>
<li>Installation</li>
</ul>
<p>After unzipping the file, follow the instructions in the README file that's included in the download package. On a Mac, moving the package to your preferred location and adding the rtg binary to your path to make it available on the command line is sufficient to complete the installation process. </p>
<ul>
<li>Testing</li>
</ul>
<p>Open a shell and run:  </p>
<pre><code class="pre_md">rtg</code class="pre_md"></pre>
<p>This should print out some usage information as well as a list of commands. As stated, to use the RTG tools you will always build your command lines like this: </p>
<pre><code class="pre_md">rtg &lt;command&gt; [options] </code class="pre_md"></pre>
<p>This means you first make the call to the binary (<code>rtg</code>), then you specify which command (method) you wish to use (e.g. <code>vcfeval</code>) then any options (<em>i.e.</em> arguments such as input files or parameters) used by the program to perform that command. This is a similar convention as used by BWA.</p>
<p>We will use RTG Tools’s modules <code>vcfeval</code> and <code>rocplot</code>. You'll find a PDF file named RTGOperationsManual.pdf containing detailed documentation included in the download packge. For our workshops, the relevant pages are pages 38–42 (for <code>vcfeval</code>) and pages 44–46 (for <code>rocplot</code>). </p>
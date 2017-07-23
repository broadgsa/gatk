## (howto) Install all software packages required to follow the GATK Best Practices.

http://gatkforums.broadinstitute.org/gatk/discussion/2899/howto-install-all-software-packages-required-to-follow-the-gatk-best-practices

<h4>Objective</h4>
<p>Install all software packages required to follow the GATK Best Practices. </p>
<h4>Prerequisites</h4>
<p>To follow these instructions, you will need to have a basic understanding of the meaning of the following words and command-line operations. If you are unfamiliar with any of the following, you should consult a more experienced colleague or your systems administrator if you have one. There are also many good online tutorials you can use to learn the necessary notions.</p>
<ul>
<li>Basic Unix environment commands </li>
<li>Binary / Executable </li>
<li>Compiling a binary </li>
<li>Adding a binary to your path </li>
<li>Command-line shell, terminal or console </li>
<li>Software library</li>
</ul>
<p>You will also need to have access to an ANSI compliant C++ compiler and the tools needed for normal compilations (make, shell, the standard library, tar, gunzip). These tools are usually pre-installed on Linux/Unix systems. <strong>On MacOS X, you may need to install the MacOS Xcode tools.</strong> See <a href="https://developer.apple.com/xcode/">https://developer.apple.com/xcode/</a> for relevant information and software downloads. The XCode tools are free but an AppleID may be required to download them.</p>
<p>Starting with version 3.6, the GATK requires Java Runtime Environment version 1.8 (Java 8). Previous versions down to 2.6 required JRE 1.7, and earlier versions required 1.6. All Linux/Unix and MacOS X systems should have a JRE pre-installed, but the version may vary. To test your Java version, run the following command in the shell: </p>
<pre><code class="pre_md">java -version </code class="pre_md"></pre>
<p>This should return a message along the lines of ”java version 1.8.0_25” as well as some details on the Runtime Environment (JRE) and Virtual Machine (VM). If you have a version that does not match the requirements stated above for the version of GATK you are running, the GATK may not run correctly or at all. The simplest solution is to install an additional JRE and specify which you want to use at the command-line. To find out how to do so, you should seek help from your systems administrator. </p>
<h4>Software packages</h4>
<ol>
<li>BWA</li>
<li>SAMtools</li>
<li>Picard</li>
<li>Genome Analysis Toolkit (GATK) </li>
<li>IGV  </li>
<li>RStudio IDE and R libraries ggplot2 and gsalib  </li>
</ol>
<p><em>Note that the version numbers of packages you download may be different than shown in the instructions below. If so, please adapt the number accordingly in the commands.</em></p>
<hr />
<h3>1. BWA</h3>
<p>Read the overview of the BWA software on the <a href="http://bio-bwa.sourceforge.net/">BWA project homepage</a>, then download the <a href="http://sourceforge.net/projects/bio-bwa/files/">latest version of the software package</a>. </p>
<ul>
<li>Installation</li>
</ul>
<p>Unpack the tar file using:   </p>
<pre><code class="pre_md">tar xvzf bwa-0.7.12.tar.bz2 </code class="pre_md"></pre>
<p>This will produce a directory called <code>bwa-0.7.12</code> containing the files necessary to compile the BWA binary. Move to this directory and compile using:   </p>
<pre><code class="pre_md">cd bwa-0.7.12
make</code class="pre_md"></pre>
<p>The compiled binary is called <code>bwa</code>. You should find it within the same folder (<code>bwa-0.7.12</code> in this example). You may also find other compiled binaries; at time of writing, a second binary called <code>bwamem-lite</code> is also included. You can disregard this file for now. Finally, just add the BWA binary to your path to make it available on the command line. This completes the installation process. </p>
<ul>
<li>Testing</li>
</ul>
<p>Open a shell and run: </p>
<pre><code class="pre_md">bwa </code class="pre_md"></pre>
<p>This should print out some version and author information as well as a list of commands. As the <strong>Usage</strong> line states, to use BWA you will always build your command lines like this: </p>
<pre><code class="pre_md">bwa &lt;command&gt; [options] </code class="pre_md"></pre>
<p>This means you first make the call to the binary (<code>bwa</code>), then you specify which command (method) you wish to use (e.g. <code>index</code>) then any options (<em>i.e.</em> arguments such as input files or parameters) used by the program to perform that command.</p>
<hr />
<h3>2. SAMtools</h3>
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
<h3>3. Picard</h3>
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
<h3>4. Genome Analysis Toolkit (GATK)</h3>
<p>Hopefully if you're reading this, you're already acquainted with the <a href="http://www.broadinstitute.org/gatk/about">purpose of the GATK</a>, so go ahead and download the <a href="http://www.broadinstitute.org/gatk/download">latest version of the software package</a>. </p>
<p>In order to access the downloads, you need to register for a free account on the <a href="http://gatkforums.broadinstitute.org/">GATK support forum</a>. You will also need to read and accept the license agreement before downloading the GATK software package. Note that if you intend to use the GATK for commercial purposes, you will need to purchase a license. See the <a href="https://www.broadinstitute.org/gatk/about/#licensing">licensing page</a> for an overview of the commercial licensing conditions. </p>
<ul>
<li>Installation</li>
</ul>
<p>Unpack the tar file using:  </p>
<pre><code class="pre_md">tar xjf GenomeAnalysisTK-3.3-0.tar.bz2 </code class="pre_md"></pre>
<p>This will produce a directory called <code>GenomeAnalysisTK-3.3-0</code> containing the GATK jar file, which is called <code>GenomeAnalysisTK.jar</code>, as well as a directory of example files called <code>resources</code>. GATK tools are distributed as a single pre-compiled Java executable so there is no need to compile them. Just like we discussed for Picard, it's not possible to add the GATK to your path, but you can set up a shortcut to the jar file using environment variables as described above. </p>
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
<h3>5. IGV</h3>
<p>The Integrated Genomics Viewer is a genome browser that allows you to view BAM, VCF and other genomic file information in context. It has a graphical user interface that is very easy to use, and can be downloaded for free (though registration is required) from <a href="https://www.broadinstitute.org/igv/home">this website</a>. We encourage you to read through IGV's very helpful <a href="https://www.broadinstitute.org/software/igv/UserGuide">user guide</a>, which includes many detailed tutorials that will help you use the program most effectively. </p>
<hr />
<h3>6. RStudio IDE and R libraries ggplot2 and gsalib</h3>
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
<p><strong>Important note</strong></p>
<p>If you are using a recent version of <code>ggplot2</code> and a version of GATK older than 3.2, you may encounter an error when trying to generate the BQSR or VQSR recalibration plots. This is because until recently our scripts were still using an older version of certain <code>ggplot2</code> functions. This has been fixed in GATK 3.2, so you should either upgrade your version of GATK (recommended) or downgrade your version of ggplot2.  If you experience further issues generating the BQSR recalibration plots, please see <a href="http://www.broadinstitute.org/gatk/guide/article?id=4294">this tutorial</a>. </p>
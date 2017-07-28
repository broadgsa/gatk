## What are the prerequisites for running GATK?

http://gatkforums.broadinstitute.org/gatk/discussion/1852/what-are-the-prerequisites-for-running-gatk

<h3>1. Operating system</h3>
<p>The GATK runs natively on most if not all flavors of UNIX, which includes MacOSX, Linux and BSD. It is possible to get it running on Windows using Cygwin, but we don't provide any support nor instructions for that.</p>
<h3>2. Java 7 / 1.7</h3>
<p>The GATK is a Java-based program, so you'll need to have Java installed on your machine. The Java version should be at 1.7 (at this time we don't officially support 1.8, and 1.6 no longer works). You can check what version you have by typing <code>java -version</code> at the command line. <a href="http://www.broadinstitute.org/gatk/guide/article?id=1200">This article</a> has some more details about what to do if you don't have the right version. Note that at this time we only support the Sun/Oracle Java JDK; OpenJDK is not supported. </p>
<h3>4. R dependencies</h3>
<p>Some of the GATK tools produce plots using R, so if you want to get the plots you'll need to have R and Rscript installed, as well as several R libraries. Full details can be found in the <a href="http://www.broadinstitute.org/gatk/guide/article?id=2899">Tutorial on installing required software</a>.</p>
<h3>3. Familiarity with command-line programs</h3>
<p>The GATK does not have a Graphical User Interface (GUI). You don't open it by clicking on the <code>.jar</code> file; you have to use the Console (or Terminal) to input commands. If this is all new to you, we recommend you first learn about that and follow some <a href="http://lifehacker.com/5633909/who-needs-a-mouse-learn-to-use-the-command-line-for-almost-anything">online tutorials</a> before trying to use the GATK. It's not difficult but you'll need to learn some jargon and get used to living without a mouse. Trust us, it's a liberating experience :)</p>
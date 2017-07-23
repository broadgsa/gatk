## What is "Phone Home" and how does it affect me?

http://gatkforums.broadinstitute.org/gatk/discussion/1250/what-is-phone-home-and-how-does-it-affect-me

<p>In GATK versions produced between September 2010 and May 2016, the GATK had a &quot;Phone Home&quot; usage reporting feature that sent us information about each GATK run via the Broad filesystem (within the Broad) and Amazon's S3 cloud storage service (outside the Broad). This feature was enabled by default and required a key to be disabled (for running offline or for regulatory reasons).</p>
<p><strong>The Phone Home feature was removed in version 3.6.</strong> Keys are no longer necessary, so if you had one, you can stop using it. We do not expect that including Phone Home arguments in GATK command lines would cause any errors (so this should not break any scripts), but let us know if you run into any trouble.</p>
<p>Note that keys remain necessary for disabling Phone Home in older versions of GATK. See further below  for details on how to obtain a key. </p>
<hr />
<h3>How Phone Home helped development</h3>
<p>At the time, the information provided by the Phone Home feature was critical in driving improvements to the GATK:</p>
<ul>
<li>By recording detailed information about each error that occurs, it enabled GATK developers to <strong>identify and fix previously-unknown bugs</strong> in the GATK. </li>
<li>It allowed us to better understand how the GATK is used in practice and <strong>adjust our documentation and development goals</strong> for common use cases.</li>
<li>It gave us a picture of <strong>which versions</strong> of the GATK are in use over time, and how successful we've been at encouraging users to migrate from obsolete or broken versions of the GATK to newer, improved versions.</li>
<li>It told us <strong>which tools</strong> were most commonly used, allowing us to monitor the adoption of newly-released tools and abandonment of outdated tools.</li>
<li>It provided us with a sense of the <strong>overall size of our user base</strong> and the major organizations/institutions using the GATK.</li>
</ul>
<hr />
<h3>What information was sent to us</h3>
<p>Below are two example GATK Run Reports showing exactly what information is sent to us each time the GATK phones home.</p>
<h4>A successful run:</h4>
<pre><code class="pre_md">&lt;GATK-run-report&gt;
    &lt;id&gt;D7D31ULwTSxlAwnEOSmW6Z4PawXwMxEz&lt;/id&gt;
    &lt;start-time&gt;2012/03/10 20.21.19&lt;/start-time&gt;
    &lt;end-time&gt;2012/03/10 20.21.19&lt;/end-time&gt;
    &lt;run-time&gt;0&lt;/run-time&gt;
    &lt;walker-name&gt;CountReads&lt;/walker-name&gt;
    &lt;svn-version&gt;1.4-483-g63ecdb2&lt;/svn-version&gt;
    &lt;total-memory&gt;85000192&lt;/total-memory&gt;
    &lt;max-memory&gt;129957888&lt;/max-memory&gt;
    &lt;user-name&gt;depristo&lt;/user-name&gt;
    &lt;host-name&gt;10.0.1.10&lt;/host-name&gt;
    &lt;java&gt;Apple Inc.-1.6.0_26&lt;/java&gt;
    &lt;machine&gt;Mac OS X-x86_64&lt;/machine&gt;
    &lt;iterations&gt;105&lt;/iterations&gt;
&lt;/GATK-run-report&gt;</code class="pre_md"></pre>
<h4>A run where an exception has occurred:</h4>
<pre><code class="pre_md">&lt;GATK-run-report&gt;
   &lt;id&gt;yX3AnltsqIlXH9kAQqTWHQUd8CQ5bikz&lt;/id&gt;   
   &lt;exception&gt;
      &lt;message&gt;Failed to parse Genome Location string: 20:10,000,000-10,000,001x&lt;/message&gt;
      &lt;stacktrace class="java.util.ArrayList"&gt; 
         &lt;string&gt;org.broadinstitute.sting.utils.GenomeLocParser.parseGenomeLoc(GenomeLocParser.java:377)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.utils.interval.IntervalUtils.parseIntervalArguments(IntervalUtils.java:82)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.commandline.IntervalBinding.getIntervals(IntervalBinding.java:106)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.gatk.GenomeAnalysisEngine.loadIntervals(GenomeAnalysisEngine.java:618)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.gatk.GenomeAnalysisEngine.initializeIntervals(GenomeAnalysisEngine.java:585)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.gatk.GenomeAnalysisEngine.execute(GenomeAnalysisEngine.java:231)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.gatk.CommandLineExecutable.execute(CommandLineExecutable.java:128)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.commandline.CommandLineProgram.start(CommandLineProgram.java:236)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.commandline.CommandLineProgram.start(CommandLineProgram.java:146)&lt;/string&gt;
         &lt;string&gt;org.broadinstitute.sting.gatk.CommandLineGATK.main(CommandLineGATK.java:92)&lt;/string&gt;
      &lt;/stacktrace&gt;
      &lt;cause&gt;
         &lt;message&gt;Position: &amp;apos;10,000,001x&amp;apos; contains invalid chars.&lt;/message&gt;
         &lt;stacktrace class="java.util.ArrayList"&gt;
            &lt;string&gt;org.broadinstitute.sting.utils.GenomeLocParser.parsePosition(GenomeLocParser.java:411)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.utils.GenomeLocParser.parseGenomeLoc(GenomeLocParser.java:374)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.utils.interval.IntervalUtils.parseIntervalArguments(IntervalUtils.java:82)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.commandline.IntervalBinding.getIntervals(IntervalBinding.java:106)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.gatk.GenomeAnalysisEngine.loadIntervals(GenomeAnalysisEngine.java:618)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.gatk.GenomeAnalysisEngine.initializeIntervals(GenomeAnalysisEngine.java:585)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.gatk.GenomeAnalysisEngine.execute(GenomeAnalysisEngine.java:231)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.gatk.CommandLineExecutable.execute(CommandLineExecutable.java:128)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.commandline.CommandLineProgram.start(CommandLineProgram.java:236)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.commandline.CommandLineProgram.start(CommandLineProgram.java:146)&lt;/string&gt;
            &lt;string&gt;org.broadinstitute.sting.gatk.CommandLineGATK.main(CommandLineGATK.java:92)&lt;/string&gt;
         &lt;/stacktrace&gt;
         &lt;is-user-exception&gt;false&lt;/is-user-exception&gt;
      &lt;/cause&gt;
      &lt;is-user-exception&gt;true&lt;/is-user-exception&gt;
   &lt;/exception&gt;
   &lt;start-time&gt;2012/03/10 20.19.52&lt;/start-time&gt;
   &lt;end-time&gt;2012/03/10 20.19.52&lt;/end-time&gt;
   &lt;run-time&gt;0&lt;/run-time&gt;
   &lt;walker-name&gt;CountReads&lt;/walker-name&gt;
   &lt;svn-version&gt;1.4-483-g63ecdb2&lt;/svn-version&gt;
   &lt;total-memory&gt;85000192&lt;/total-memory&gt;
   &lt;max-memory&gt;129957888&lt;/max-memory&gt;
   &lt;user-name&gt;depristo&lt;/user-name&gt;
   &lt;host-name&gt;10.0.1.10&lt;/host-name&gt;
   &lt;java&gt;Apple Inc.-1.6.0_26&lt;/java&gt;
   &lt;machine&gt;Mac OS X-x86_64&lt;/machine&gt;
   &lt;iterations&gt;0&lt;/iterations&gt;
&lt;/GATK-run-report&gt;</code class="pre_md"></pre>
<p><strong>Note that as of GATK 1.5 we no longer collected information about the command-line executed, the working directory, or tmp directory.</strong></p>
<hr />
<h3>Disabling Phone Home</h3>
<p>Versions of GATK older than 3.6 attempted to &quot;phone home&quot; as a normal part of each run. However, we recognized that some of our users need to run the GATK with the Phone Home disabled. To enable this, we provided an option (<code>-et NO_ET</code> )  in GATK 1.5 and later to disable the Phone Home feature. To use this option, you need to contact us to request a key. Instructions for doing so are below.</p>
<h4>How to obtain and use a GATK key</h4>
<p>To obtain a GATK key, please fill out the <a href="http://www.broadinstitute.org/gatk/request-key">request form</a>. </p>
<p>Running the GATK with a key is simple: you just need to append a <code>-K your.key</code> argument to your customary command line, where <code>your.key</code> is the path to the key file you obtained from us:</p>
<pre><code class="pre_md">java -jar dist/GenomeAnalysisTK.jar \
    -T PrintReads \
    -I public/testdata/exampleBAM.bam \
    -R public/testdata/exampleFASTA.fasta \
    -et NO_ET \
    -K your.key</code class="pre_md"></pre>
<p>The <code>-K</code> argument is only necessary when running the GATK with the <code>NO_ET</code> option.</p>
<h4>Troubleshooting key-related problems</h4>
<ul>
<li>Corrupt/Unreadable/Revoked Keys</li>
</ul>
<p>If you get an error message from the GATK saying that your key is corrupt, unreadable, or has been revoked, please apply for a new key.</p>
<ul>
<li>GATK Public Key Not Found</li>
</ul>
<p>If you get an error message stating that the GATK public key could not be located or read, then something is likely wrong with your build of the GATK. If you're running the binary release, try <a href="http://www.broadinstitute.org/gatk/download">downloading</a> it again. If you're compiling from source, try re-compiling. If all else fails, please ask for help on our <a href="http://gatkforums.broadinstitute.org/">community forum</a>.</p>
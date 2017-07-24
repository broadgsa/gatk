## Writing unit tests for walkers

http://gatkforums.broadinstitute.org/gatk/discussion/1339/writing-unit-tests-for-walkers

<h2>1. Testing core walkers is critical</h2>
<p>Most GATK walkers are really too complex to easily test using the standard unit test framework.  It's just not feasible to make artificial read piles and then extrapolate from simple tests passing whether the system as a whole is working correctly.  However, we need some way to determine whether changes to the core of the GATK are altering the expected output of complex walkers like BaseRecalibrator or SingleSampleGenotyper.   In additional to correctness, we want to make sure that the performance of key walkers isn't degrading over time, so that calling snps, cleaning indels, etc., isn't slowly creeping down over time.  Since we are now using a bamboo server to automatically build and run unit tests (as well as measure their runtimes) we want to put as many good walker tests into the test framework so we capture performance metrics over time.</p>
<h2>2. The WalkerTest framework</h2>
<p>To make this testing process easier, we've created a <code>WalkerTest</code> framework that lets you invoke the GATK using command-line GATK commands in the <code>JUnit</code> system and test for changes in your output files by comparing the current ant build results to previous run via an MD5 sum.  It's a bit coarse grain, but it will work to ensure that changes to key walkers are detected quickly by the system, and authors can either update the expected MD5s or go track down bugs.</p>
<p>The system is fairly straightforward to use.  Ultimately we will end up with <code>JUnit</code> style tests in the unit testing structure.  In the piece of code below, we have a piece of code that checks the MD5 of the SingleSampleGenotyper's GELI text output at LOD 3 and LOD 10.  </p>
<pre><code class="pre_md">package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

public class SingleSampleGenotyperTest extends WalkerTest {
    @Test
    public void testLOD() {
        HashMap&lt;Double, String&gt; e = new HashMap&lt;Double, String&gt;();
        e.put( 10.0, "e4c51dca6f1fa999f4399b7412829534" );
        e.put( 3.0, "d804c24d49669235e3660e92e664ba1a" );

        for ( Map.Entry&lt;Double, String&gt; entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                   "-T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s --variant_output_format GELI -L 1:10,000,000-11,000,000 -m EMPIRICAL -lod " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest("testLOD", spec);
        }
    }
}</code class="pre_md"></pre>
<p>The fundamental piece here is to inherit from <code>WalkerTest</code>.  This gives you access to the <code>executeTest()</code> function that consumes a <code>WalkerTestSpec</code>:</p>
<pre><code class="pre_md">    public WalkerTestSpec(String args, int nOutputFiles, List&lt;String&gt; md5s)</code class="pre_md"></pre>
<p>The <code>WalkerTestSpec</code> takes regular, command-line style GATK arguments describing what you want to run, the number of output files the walker will generate, and your expected MD5s for each of these output files.  The args string can contain <code>%s String.format</code> specifications, and for each of the <code>nOutputFiles</code>, the <code>executeTest()</code> function will (1) generate a <code>tmp</code> file for output and (2) call <code>String.format</code> on your args to fill in the tmp output files in your arguments string.  For example, in the above argument string <code>varout</code> is followed by <code>%s</code>, so our single SingleSampleGenotyper output is the variant output file.</p>
<h2>3. Example output</h2>
<p>When you add a <code>walkerTest</code> inherited unit test to the GATK, and then <code>build test</code>, you'll see output that looks like:</p>
<pre><code class="pre_md">[junit] WARN  13:29:50,068 WalkerTest - -------------------------------------------------------------------------------- 
[junit] WARN  13:29:50,068 WalkerTest - -------------------------------------------------------------------------------- 
[junit] WARN  13:29:50,069 WalkerTest - Executing test testLOD with GATK arguments: -T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout /tmp/walktest.tmp_param.05524470250256847817.tmp --variant_output_format GELI -L 1:10,000,000-11,000,000 -m EMPIRICAL -lod 3.0
[junit]  
[junit] WARN  13:29:50,069 WalkerTest - Executing test testLOD with GATK arguments: -T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout /tmp/walktest.tmp_param.05524470250256847817.tmp --variant_output_format GELI -L 1:10,000,000-11,000,000 -m EMPIRICAL -lod 3.0
[junit]  
[junit] WARN  13:30:39,407 WalkerTest - Checking MD5 for /tmp/walktest.tmp_param.05524470250256847817.tmp [calculated=d804c24d49669235e3660e92e664ba1a, expected=d804c24d49669235e3660e92e664ba1a] 
[junit] WARN  13:30:39,407 WalkerTest - Checking MD5 for /tmp/walktest.tmp_param.05524470250256847817.tmp [calculated=d804c24d49669235e3660e92e664ba1a, expected=d804c24d49669235e3660e92e664ba1a] 
[junit] WARN  13:30:39,408 WalkerTest -   =&gt; testLOD PASSED 
[junit] WARN  13:30:39,408 WalkerTest -   =&gt; testLOD PASSED 
[junit] WARN  13:30:39,409 WalkerTest - -------------------------------------------------------------------------------- 
[junit] WARN  13:30:39,409 WalkerTest - -------------------------------------------------------------------------------- 
[junit] WARN  13:30:39,409 WalkerTest - Executing test testLOD with GATK arguments: -T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout /tmp/walktest.tmp_param.03852477489430798188.tmp --variant_output_format GELI -L 1:10,000,000-11,000,000 -m EMPIRICAL -lod 10.0
[junit]  
[junit] WARN  13:30:39,409 WalkerTest - Executing test testLOD with GATK arguments: -T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout /tmp/walktest.tmp_param.03852477489430798188.tmp --variant_output_format GELI -L 1:10,000,000-11,000,000 -m EMPIRICAL -lod 10.0
[junit]  
[junit] WARN  13:31:30,213 WalkerTest - Checking MD5 for /tmp/walktest.tmp_param.03852477489430798188.tmp [calculated=e4c51dca6f1fa999f4399b7412829534, expected=e4c51dca6f1fa999f4399b7412829534] 
[junit] WARN  13:31:30,213 WalkerTest - Checking MD5 for /tmp/walktest.tmp_param.03852477489430798188.tmp [calculated=e4c51dca6f1fa999f4399b7412829534, expected=e4c51dca6f1fa999f4399b7412829534] 
[junit] WARN  13:31:30,213 WalkerTest -   =&gt; testLOD PASSED 
[junit] WARN  13:31:30,213 WalkerTest -   =&gt; testLOD PASSED 
[junit] WARN  13:31:30,214 SingleSampleGenotyperTest -  
[junit] WARN  13:31:30,214 SingleSampleGenotyperTest -  </code class="pre_md"></pre>
<h2>4. Recommended location for GATK testing data</h2>
<p>We keep all of the permenant GATK testing data in:</p>
<pre><code class="pre_md">/humgen/gsa-scr1/GATK_Data/Validation_Data/</code class="pre_md"></pre>
<p>A good set of data to use for walker testing is the CEU daughter data from 1000 Genomes:</p>
<pre><code class="pre_md">gsa2 ~/dev/GenomeAnalysisTK/trunk &gt; ls -ltr /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_1*.bam /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_1*.calls
-rw-rw-r--+ 1 depristo wga  51M 2009-09-03 07:56 /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam
-rw-rw-r--+ 1 depristo wga 185K 2009-09-04 13:21 /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls
-rw-rw-r--+ 1 depristo wga 164M 2009-09-04 13:22 /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls
-rw-rw-r--+ 1 depristo wga  24M 2009-09-04 15:00 /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam
-rw-rw-r--+ 1 depristo wga  12M 2009-09-04 15:01 /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.454.bam
-rw-r--r--+ 1 depristo wga  91M 2009-09-04 15:02 /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam</code class="pre_md"></pre>
<h2>5. Test dependencies</h2>
<p>The tests depend on a variety of input files, that are generally constrained to three mount points on the internal Broad network:</p>
<pre><code class="pre_md">*/seq/
*/humgen/1kg/
*/humgen/gsa-hpprojects/GATK/Data/Validation_Data/</code class="pre_md"></pre>
<p>To run the unit and integration tests you'll have to have access to these files.  They may have different mount points on your machine (say, if you're running remotely over the VPN and have mounted the directories on your own machine).</p>
<h2>6. MD5 database and comparing MD5 results</h2>
<p>Every file that generates an MD5 sum as part of the WalkerTest framework will be copied to <code>&lt;MD5&gt;. integrationtest</code> in the <code>integrationtests</code> subdirectory of the GATK trunk.  This MD5 database of results enables you to easily examine the results of an integration test as well as compare the results of a test before/after a code change.  For example, below is an example test for the UnifiedGenotyper that, due to a code change, where the output VCF differs from the VCF with the expected MD5 value in the test code itself.  The test provides provides the path to the two results files as well as a diff command to compare expected to the observed MD5:</p>
<pre><code class="pre_md">[junit] --------------------------------------------------------------------------------    
[junit] Executing test testParameter[-genotype] with GATK arguments: -T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout /tmp/walktest.tmp_param.05997727998894311741.tmp -L 1:10,000,000-10,010,000 -genotype    
[junit] ##### MD5 file is up to date: integrationtests/ab20d4953b13c3fc3060d12c7c6fe29d.integrationtest    
[junit] Checking MD5 for /tmp/walktest.tmp_param.05997727998894311741.tmp [calculated=ab20d4953b13c3fc3060d12c7c6fe29d, expected=0ac7ab893a3f550cb1b8c34f28baedf6]    
[junit] ##### Test testParameter[-genotype] is going fail #####    
[junit] ##### Path to expected   file (MD5=0ac7ab893a3f550cb1b8c34f28baedf6): integrationtests/0ac7ab893a3f550cb1b8c34f28baedf6.integrationtest    
[junit] ##### Path to calculated file (MD5=ab20d4953b13c3fc3060d12c7c6fe29d): integrationtests/ab20d4953b13c3fc3060d12c7c6fe29d.integrationtest    
[junit] ##### Diff command: diff integrationtests/0ac7ab893a3f550cb1b8c34f28baedf6.integrationtest integrationtests/ab20d4953b13c3fc3060d12c7c6fe29d.integrationtest</code class="pre_md"></pre>
<p>Examining the diff we see a few lines that have changed the <code>DP</code> count in the new code</p>
<pre><code class="pre_md">&gt; diff integrationtests/0ac7ab893a3f550cb1b8c34f28baedf6.integrationtest integrationtests/ab20d4953b13c3fc3060d12c7c6fe29d.integrationtest  | head
385,387c385,387
&lt; 1     10000345        .       A       .       106.54  .       AN=2;DP=33;Dels=0.00;MQ=89.17;MQ0=0;SB=-10.00   GT:DP:GL:GQ     0/0:25:-0.09,-7.57,-75.74:74.78
&lt; 1     10000346        .       A       .       103.75  .       AN=2;DP=31;Dels=0.00;MQ=88.85;MQ0=0;SB=-10.00   GT:DP:GL:GQ     0/0:24:-0.07,-7.27,-76.00:71.99
&lt; 1     10000347        .       A       .       109.79  .       AN=2;DP=31;Dels=0.00;MQ=88.85;MQ0=0;SB=-10.00   GT:DP:GL:GQ     0/0:26:-0.05,-7.85,-84.74:78.04
---
&gt; 1     10000345        .       A       .       106.54  .       AN=2;DP=32;Dels=0.00;MQ=89.50;MQ0=0;SB=-10.00   GT:DP:GL:GQ     0/0:25:-0.09,-7.57,-75.74:74.78
&gt; 1     10000346        .       A       .       103.75  .       AN=2;DP=30;Dels=0.00;MQ=89.18;MQ0=0;SB=-10.00   GT:DP:GL:GQ     0/0:24:-0.07,-7.27,-76.00:71.99
&gt; 1     10000347        .       A       .       109.79  .       AN=2;DP=30;Dels=0.00;MQ=89.18;MQ0=0;SB=-10.00   GT:DP:GL:GQ     0/0:26:-0.05,-7.85,-84.74:78</code class="pre_md"></pre>
<p>Whether this is the expected change is up to you to decide, but the system makes it as easy as possible to see the consequences of your code change.</p>
<h2>7. Testing for Exceptions</h2>
<p>The walker test framework supports an additional syntax for ensuring that a particular java Exception is thrown when a walker executes using a simple alternate version of the <code>WalkerSpec</code> object.  Rather than specifying the MD5 of the result, you can provide a single subclass of <code>Exception.class</code> and the testing framework will ensure that when the walker runs an instance (class or subclass) of your expected exception is thrown.  The system also flags if no exception is thrown.</p>
<p>For example, the following code tests that the GATK can detect and error out when incompatible VCF and FASTA files are given:</p>
<pre><code class="pre_md">@Test public void fail8() { executeTest("hg18lex-v-b36", test(lexHG18, callsB36)); }

private WalkerTest.WalkerTestSpec test(String ref, String vcf) {
    return new WalkerTest.WalkerTestSpec("-T VariantsToTable -M 10 -B:two,vcf "
            + vcf + " -F POS,CHROM -R "
            + ref +  " -o %s",
            1, UserException.IncompatibleSequenceDictionaries.class);

}</code class="pre_md"></pre>
<p>During the integration test this looks like:</p>
<pre><code class="pre_md">[junit] Executing test hg18lex-v-b36 with GATK arguments: -T VariantsToTable -M 10 -B:two,vcf /humgen/gsa-hpprojects/GATK/data/Validation_Data/lowpass.N3.chr1.raw.vcf -F POS,CHROM -R /humgen/gsa-hpprojects/GATK/data/Validation_Data/lexFasta/lex.hg18.fasta -o /tmp/walktest.tmp_param.05541601616101756852.tmp -l WARN -et NO_ET
[junit]    [junit] Wanted exception class org.broadinstitute.sting.utils.exceptions.UserException$IncompatibleSequenceDictionaries, saw class org.broadinstitute.sting.utils.exceptions.UserException$IncompatibleSequenceDictionaries
[junit]   =&gt; hg18lex-v-b36 PASSED</code class="pre_md"></pre>
<h2>8. Miscellaneous information</h2>
<ul>
<li>
<p>Please do not put any extremely long tests in the regular <code>ant build test</code> target.  We are currently splitting the system into fast and slow tests so that unit tests can be run in \&lt; 3 minutes while saving a test target for long-running regression tests.  More information on that will be posted. </p>
</li>
<li>
<p>An expected MG5 string of <code>""</code> means don't check for equality between the calculated and expected MD5s.  Useful if you are just writing a new test and don't know the true output.</p>
</li>
<li>
<p>Overload <code>parameterize() { return true; }</code> if you want the system to just run your calculations, not throw an error if your MD5s don't match, across all tests</p>
</li>
<li>
<p>If your tests all of a sudden stop giving equality MD5s, you can just (1) look at the <code>.tmp</code> output files directly or (2) grab the printed GATK command-line options and explore what is happening.</p>
</li>
<li>
<p>You can always run a GATK walker on the command line and then run md5sum on its output files to obtain, outside of the testing framework, the MD5 expected results.</p>
</li>
<li>Don't worry about the duplication of lines in the output ; it's just an annoyance of having two global loggers. Eventually we'll bug fix this away.</li>
</ul>
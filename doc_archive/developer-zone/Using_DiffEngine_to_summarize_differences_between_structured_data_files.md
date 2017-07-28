## Using DiffEngine to summarize differences between structured data files

http://gatkforums.broadinstitute.org/gatk/discussion/1299/using-diffengine-to-summarize-differences-between-structured-data-files

<h3>1. What is DiffEngine?</h3>
<p>DiffEngine is a summarizing difference engine that allows you to compare two structured files -- such as BAMs and VCFs -- to find what are the differences between them.  This is primarily useful in regression testing or optimization, where you want to ensure that the differences are those that you expect and not any others.  </p>
<h3>2. The summarized differences</h3>
<p>The GATK contains a summarizing difference engine called DiffEngine that compares hierarchical data structures to emit:</p>
<ul>
<li>
<p>A list of specific differences between the two data structures.  This is similar to saying the value in field A in record 1 in file F differs from the value in field A in record 1 in file G.  </p>
</li>
<li>A summarized list of differences ordered by frequency of the difference.  This output is similar to saying field A differed in 50 records between files F and G.</li>
</ul>
<h3>3. The DiffObjects walker</h3>
<p>The GATK contains a private walker called DiffObjects that allows you access to the DiffEngine capabilities on the command line.  Simply provide the walker with the master and test files and it will emit summarized differences for you.</p>
<h3>4. Understanding the output</h3>
<p>The DiffEngine system compares to two hierarchical data structures for specific differences in the values of named nodes.  Suppose I have two trees:</p>
<pre><code class="pre_md">Tree1=(A=1 B=(C=2 D=3)) 
Tree2=(A=1 B=(C=3 D=3 E=4))
Tree3=(A=1 B=(C=4 D=3 E=4))</code class="pre_md"></pre>
<p>where every node in the tree is named, or is a raw value (here all leaf values are integers).  The DiffEngine traverses these data structures by name, identifies equivalent nodes by fully qualified names (<code>Tree1.A</code> is distinct from <code>Tree2.A</code>, and determines where their values are equal (<code>Tree1.A=1</code>, <code>Tree2.A=1</code>, so they are).  </p>
<p>These itemized differences are listed as:</p>
<pre><code class="pre_md">Tree1.B.C=2 != Tree2.B.C=3
Tree1.B.C=2 != Tree3.B.C=4
Tree2.B.C=3 != Tree3.B.C=4
Tree1.B.E=MISSING != Tree2.B.E=4</code class="pre_md"></pre>
<p>This conceptually very similar to the output of the unix command line tool <code>diff</code>.  What's nice about DiffEngine though is that it computes similarity among the itemized differences and displays the count of differences names in the system.  In the above example, the field <code>C</code> is not equal three times, while the missing <code>E</code> in <code>Tree1</code> occurs only once.  So the summary is:</p>
<pre><code class="pre_md">*.B.C : 3
*.B.E : 1</code class="pre_md"></pre>
<p>where the <code>*</code> operator indicates that any named field matches.  This output is sorted by counts, and provides an immediate picture of the commonly occurring differences between the files.  </p>
<p>Below is a detailed example of two VCF fields that differ because of a bug in the <code>AC</code>, <code>AF</code>, and <code>AN</code> counting routines, detected by the <code>integrationtest</code> integration (more below).  You can see that in the although there are many specific instances of these differences between the two files, the summarized differences provide an immediate picture that the <code>AC</code>, <code>AF</code>, and <code>AN</code> fields are the major causes of the differences.</p>
<pre><code class="pre_md">[testng] path                                                              count
[testng] *.*.*.AC                                                         6
[testng] *.*.*.AF                                                         6
[testng] *.*.*.AN                                                         6
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000000.AC  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000000.AF  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000000.AN  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000117.AC  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000117.AF  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000117.AN  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000211.AC  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000211.AF  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000211.AN  1
[testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000598.AC  1</code class="pre_md"></pre>
<h3>5. Integration tests</h3>
<p>The DiffEngine codebase that supports these calculations is integrated into the <code>integrationtest</code> framework, so that when a test fails the system automatically summarizes the differences between the master MD5 file and the failing MD5 file, if it is an understood type.  When failing you will see in the integration test logs not only the basic information, but the detailed DiffEngine output.  </p>
<p>For example, in the output below I broke the GATK BAQ calculation and the integration test DiffEngine clearly identifies that all of the records differ in their <code>BQ</code> tag value in the two BAM files:</p>
<pre><code class="pre_md">/humgen/1kg/reference/human_b36_both.fasta -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam -o /var/folders/Us/UsMJ3xRrFVyuDXWkUos1xkC43FQ/-Tmp-/walktest.tmp_param.05785205687740257584.tmp -L 1:10,000,000-10,100,000 -baq RECALCULATE -et NO_ET
   [testng] WARN  22:59:22,875 TextFormattingUtils - Unable to load help text.  Help output will be sparse.
   [testng] WARN  22:59:22,875 TextFormattingUtils - Unable to load help text.  Help output will be sparse.
   [testng] ##### MD5 file is up to date: integrationtests/e5147656858fc4a5f470177b94b1fc1b.integrationtest
   [testng] Checking MD5 for /var/folders/Us/UsMJ3xRrFVyuDXWkUos1xkC43FQ/-Tmp-/walktest.tmp_param.05785205687740257584.tmp [calculated=e5147656858fc4a5f470177b94b1fc1b, expected=4ac691bde1ba1301a59857694fda6ae2]
   [testng] ##### Test testPrintReadsRecalBAQ is going fail #####
   [testng] ##### Path to expected   file (MD5=4ac691bde1ba1301a59857694fda6ae2): integrationtests/4ac691bde1ba1301a59857694fda6ae2.integrationtest
   [testng] ##### Path to calculated file (MD5=e5147656858fc4a5f470177b94b1fc1b): integrationtests/e5147656858fc4a5f470177b94b1fc1b.integrationtest
   [testng] ##### Diff command: diff integrationtests/4ac691bde1ba1301a59857694fda6ae2.integrationtest integrationtests/e5147656858fc4a5f470177b94b1fc1b.integrationtest
   [testng] ##:GATKReport.v0.1 diffences : Summarized differences between the master and test files.
   [testng] See http://www.broadinstitute.org/gsa/wiki/index.php/DiffObjectsWalker_and_SummarizedDifferences for more information
   [testng] Difference                                                                               NumberOfOccurrences
   [testng] *.*.*.BQ                                                                                 895
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:2:266:272:361.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:5:245:474:254.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:5:255:178:160.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:6:158:682:495.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:6:195:591:884.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:7:165:236:848.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:7:191:223:910.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAE_0002_FC205W7AAXX:7:286:279:434.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAF_0002_FC205Y7AAXX:2:106:516:354.BQ  1
   [testng] 4ac691bde1ba1301a59857694fda6ae2.integrationtest.-XAF_0002_FC205Y7AAXX:3:102:580:518.BQ  1
   [testng]
   [testng] Note that the above list is not comprehensive.  At most 20 lines of output, and 10 specific differences will be listed.  Please use -T DiffObjects -R public/testdata/exampleFASTA.fasta -m integrationtests/4ac691bde1ba1301a59857694fda6ae2.integrationtest -t integrationtests/e5147656858fc4a5f470177b94b1fc1b.integrationtest to explore the differences more freely</code class="pre_md"></pre>
<h3>6. Adding your own DiffableObjects to the system</h3>
<p>The system dynamically finds all classes that implement the following simple interface:</p>
<pre><code class="pre_md">public interface DiffableReader {
    @Ensures("result != null")
    /**
     * Return the name of this DiffableReader type.  For example, the VCF reader returns 'VCF' and the
     * bam reader 'BAM'
     */
    public String getName();

    @Ensures("result != null")
    @Requires("file != null")
    /**
     * Read up to maxElementsToRead DiffElements from file, and return them.
     */
    public DiffElement readFromFile(File file, int maxElementsToRead);

    /**
     * Return true if the file can be read into DiffElement objects with this reader. This should
     * be uniquely true/false for all readers, as the system will use the first reader that can read the
     * file.  This routine should never throw an exception.  The VCF reader, for example, looks at the
     * first line of the file for the ##format=VCF4.1 header, and the BAM reader for the BAM_MAGIC value
     * @param file
     * @return
     */
    @Requires("file != null")
    public boolean canRead(File file);</code class="pre_md"></pre>
<p>See the VCF and BAMDiffableReaders for example implementations.  If you extend this to a new object types both the DiffObjects walker and the <code>integrationtest</code> framework will automatically work with your new file type.</p>
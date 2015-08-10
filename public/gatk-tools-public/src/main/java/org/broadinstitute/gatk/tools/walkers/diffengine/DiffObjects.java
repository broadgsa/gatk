/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.diffengine;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.diffengine.DiffElement;
import org.broadinstitute.gatk.utils.diffengine.DiffEngine;
import org.broadinstitute.gatk.utils.diffengine.Difference;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

/**
 * A generic engine for comparing tree-structured objects
 *
 * <p>
 *      This tool compares two record-oriented files, itemizing specific difference between equivalent
 *      records in the two files.  Reports both itemized and summarized differences.
 * </p>
 *
 * <h3>What are the summarized differences and the DiffObjectsWalker?</h3>
 *
 * <p>
 *     The GATK contains a summarizing difference engine that compares hierarchical data structures to emit:
 *      <ul>
 *          <li>A list of specific differences between the two data structures.  This is similar to saying the value in field A in record 1 in file F differences from the value in field A in record 1 in file G.</li>
 *          <li>A summarized list of differences ordered by frequency of the difference.  This output is similar to saying field A in 50 records in files F and G differed.</li>
 *      </ul>
 * </p>
 *
 * <p>
 *      The GATK contains a private walker DiffObjects that allows you access to the DiffEngine capabilities on the command line.  Simply provide the walker with the master and test files and it will emit summarized differences for you.
 * </p>
 *
 * <h3>Why?</h3>
 *
 * <p>
 *      The reason for this system is that it allows you to compare two structured files -- such as BAMs and VCFs -- for common differences among them.  This is primarily useful in regression testing or optimization, where you want to ensure that the differences are those that you expect and not any others.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 *      The DiffObjectsWalker works with BAM or VCF files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 *      The DiffEngine system compares to two hierarchical data structures for specific differences in the values of named
 *      nodes.  Suppose I have two trees:
 * <pre>
 *     Tree1=(A=1 B=(C=2 D=3))
 *     Tree2=(A=1 B=(C=3 D=3 E=4))
 *     Tree3=(A=1 B=(C=4 D=3 E=4))
 * </pre>
 * <p>
 *     where every node in the tree is named, or is a raw value (here all leaf values are integers).  The DiffEngine
 *      traverses these data structures by name, identifies equivalent nodes by fully qualified names
 *      (Tree1.A is distinct from Tree2.A, and determines where their values are equal (Tree1.A=1, Tree2.A=1, so they are).
 *      These itemized differences are listed as:
 * <pre>
 *     Tree1.B.C=2 != Tree2.B.C=3
 *     Tree1.B.C=2 != Tree3.B.C=4
 *     Tree2.B.C=3 != Tree3.B.C=4
 *     Tree1.B.E=MISSING != Tree2.B.E=4
 * </pre>
 *
 * <p>
 *      This conceptually very similar to the output of the unix command line tool diff.  What's nice about DiffEngine though
 *      is that it computes similarity among the itemized differences and displays the count of differences names
 *      in the system.  In the above example, the field C is not equal three times, while the missing E in Tree1 occurs
 *      only once.  So the summary is:
 *
 * <pre>
 *     *.B.C : 3
 *     *.B.E : 1
 * </pre>
 *
 * <p>
 *      where the * operator indicates that any named field matches.  This output is sorted by counts, and provides an
 *      immediate picture of the commonly occurring differences among the files.
 * <p>
 *      Below is a detailed example of two VCF fields that differ because of a bug in the AC, AF, and AN counting routines,
 *      detected by the integrationtest integration (more below).  You can see that in the although there are many specific
 *      instances of these differences between the two files, the summarized differences provide an immediate picture that
 *      the AC, AF, and AN fields are the major causes of the differences.
 * <p>
 *
 * <pre>
 [testng] path                                                             count
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
 [testng] 64b991fd3850f83614518f7d71f0532f.integrationtest.20:10000598.AC  1
 </pre>
 *
 * <h3>Caveat</h3>
 * <p>Because this is a walker, it requires that you pass a reference file. However the reference is not actually used, so it does not matter what you pass as reference.</p>
 *
 *
 * @author Mark DePristo
 * @since 7/4/11
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class DiffObjects extends RodWalker<Integer, Integer> {
    /**
     * Writes out a file of the DiffEngine format:
     *
     *     See http://www.broadinstitute.org/gatk/guide/article?id=1299 for details.
     */
    @Output(doc="File to which results should be written")
    protected PrintStream out;

    /**
     * The master file against which we will compare test.  This is one of the two required
     * files to do the comparison.  Conceptually master is the original file contained the expected
     * results, but this doesn't currently have an impact on the calculations, but might in the future.
     */
    @Input(fullName="master", shortName="m", doc="Master file: expected results", required=true)
    File masterFile;

    /**
     * The test file against which we will compare to the master.  This is one of the two required
     * files to do the comparison.  Conceptually test is the derived file from master, but this
     * doesn't currently have an impact on the calculations, but might in the future.
     */
    @Input(fullName="test", shortName="t", doc="Test file: new results to compare to the master file", required=true)
    File testFile;

    /**
     * The engine will read at most this number of objects from each of master and test files.  This reduces
     * the memory requirements for DiffObjects but does limit you to comparing at most this number of objects
     */
    @Argument(fullName="maxObjectsToRead", shortName="motr", doc="Max. number of objects to read from the files.  -1 [default] means unlimited", required=false)
    int MAX_OBJECTS_TO_READ = -1;

    @Argument(fullName="maxRawDiffsToSummarize", shortName="maxRawDiffsToSummarize", doc="Max. number of differences to include in the summary.  -1 [default] means unlimited", required=false)
    int maxRawDiffsToSummary = -1;

    @Argument(fullName="doPairwise", shortName="doPairwise", doc="If provided, we will compute the minimum pairwise differences to summary, which can be extremely expensive", required=false)
    boolean doPairwise = false;

    /**
     * The max number of differences to display when summarizing.  For example, if there are 10M differences, but
     * maxDiffs is 10, then the comparison aborts after first ten summarized differences are shown.  Note that
     * the system shows differences sorted by frequency, so these 10 would be the most common between the two files.
     * A value of 0 means show all possible differences.
     */
    @Argument(fullName="maxDiffs", shortName="M", doc="Max. number of diffs to process", required=false)
    int MAX_DIFFS = 0;

    /**
     * The maximum number of singleton (occurs exactly once between the two files) to display when writing out
     * the summary.  Only applies if maxDiffs hasn't been exceeded.  For example, if maxDiffs is 10 and maxCount1Diffs
     * is 2 and there are 20 diffs with count > 1, then only 10 are shown, all of which have count above 1.
     */
    @Argument(fullName="maxCount1Diffs", shortName="M1", doc="Max. number of diffs occuring exactly once in the file to process", required=false)
    int MAX_COUNT1_DIFFS = 0;

    /**
     * Only differences that occur more than minCountForDiff are displayed.  For example, if minCountForDiff is 10, then
     * a difference must occur at least 10 times between the two files to be shown.
     */
    @Argument(fullName="minCountForDiff", shortName="MCFD", doc="Min number of observations for a records to display", required=false)
    int minCountForDiff = 1;

    /**
     * If provided, the system will write out the summarized, individual differences.  May lead to enormous outputs,
     * depending on how many differences are found.  Note these are not sorted in any way, so if you have 10M
     * common differences in the files, you will see 10M records, whereas the final summarize will just list the
     * difference and its count of 10M.
     */
    @Argument(fullName="showItemizedDifferences", shortName="SID", doc="Should we enumerate all differences between the files?", required=false)
    boolean showItemizedDifferences = false;

    @Argument(fullName="iterations", doc="Number of iterations to perform, should be 1 unless you are doing memory testing", required=false)
    int iterations = 1;

    DiffEngine diffEngine;

    @Override
    public void initialize() {
        this.diffEngine = new DiffEngine();
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return 0;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    @Override
    public void onTraversalDone(Integer sum) {
        if ( iterations > 1 ) {
            for ( int i = 0; i < iterations; i++ ) {
                DiffEngine.SummaryReportParams params = new DiffEngine.SummaryReportParams(out, 20, 10, 0, -1, false);
                boolean success = DiffEngine.simpleDiffFiles(masterFile, testFile, MAX_OBJECTS_TO_READ, params);
                logger.info("Iteration " + i + " success " + success);
            }
        } else {
            //out.printf("Reading master file %s%n", masterFile);
            DiffElement master = diffEngine.createDiffableFromFile(masterFile, MAX_OBJECTS_TO_READ);
            logger.info(String.format("Read %d objects", master.size()));
            //out.printf("Reading test file %s%n", testFile);
            DiffElement test = diffEngine.createDiffableFromFile(testFile, MAX_OBJECTS_TO_READ);
            logger.info(String.format("Read %d objects", test.size()));

//        out.printf("Master diff objects%n");
//        out.println(master.toString());
//        out.printf("Test diff objects%n");
//        out.println(test.toString());

            List<Difference> diffs = diffEngine.diff(master, test);
            logger.info(String.format("Done computing diff with %d differences found", diffs.size()));
            if ( showItemizedDifferences ) {
                out.printf("Itemized results%n");
                for ( Difference diff : diffs )
                    out.printf("DIFF: %s%n", diff.toString());
            }

            DiffEngine.SummaryReportParams params = new DiffEngine.SummaryReportParams(out,
                    MAX_DIFFS, MAX_COUNT1_DIFFS, minCountForDiff,
                    maxRawDiffsToSummary, doPairwise);
            params.setDescending(false);
            diffEngine.reportSummarizedDifferences(diffs, params);
            logger.info(String.format("Done summarizing differences"));
        }
    }
}

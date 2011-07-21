/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.diffengine;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

/**
 * A generic engine for comparing tree-structured objects
 *
 * Compares two record-oriented files, itemizing specific difference between equivalent
 * records in the two files.  Reports both itemized and summarized differences.
 *
 * @author Mark DePristo
 * @since 7/4/11
 * @version 0.1
 */
@Requires(value={})
public class DiffObjectsWalker extends RodWalker<Integer, Integer> {
    /**
     * Writes out a file of the DiffEngine format:
     *
     *      http://www.broadinstitute.org/gsa/wiki/index.php/DiffEngine
     */
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    /**
     * The master file against which we will compare test.  This is one of the two required
     * files to do the comparison.  Conceptually master is the original file contained the expected
     * results, but this doesn't currently have an impact on the calculations, but might in the future.
     */
    @Argument(fullName="master", shortName="m", doc="Master file: expected results", required=true)
    File masterFile;

    /**
     * The test file against which we will compare to the master.  This is one of the two required
     * files to do the comparison.  Conceptually test is the derived file from master, but this
     * doesn't currently have an impact on the calculations, but might in the future.
     */
    @Argument(fullName="test", shortName="t", doc="Test file: new results to compare to the master file", required=true)
    File testFile;

    /**
     * The engine will read at most this number of objects from each of master and test files.  This reduces
     * the memory requirements for DiffObjects but does limit you to comparing at most this number of objects
     */
    @Argument(fullName="maxObjectsToRead", shortName="motr", doc="Max. number of objects to read from the files.  -1 [default] means unlimited", required=false)
    int MAX_OBJECTS_TO_READ = -1;

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

    @Argument(fullName="testEnum", doc="X", required=false)
    TestEnum testEnum = TestEnum.ONE;

    public enum TestEnum { ONE, TWO };

    final DiffEngine diffEngine = new DiffEngine();

    @Override
    public void initialize() {

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
        out.printf("Reading master file %s%n", masterFile);
        DiffElement master = diffEngine.createDiffableFromFile(masterFile, MAX_OBJECTS_TO_READ);
        out.printf("  Read %d objects%n", master.size());
        out.printf("Reading test file %s%n", testFile);
        DiffElement test = diffEngine.createDiffableFromFile(testFile, MAX_OBJECTS_TO_READ);
        out.printf("  Read %d objects%n", test.size());

//        out.printf("Master diff objects%n");
//        out.println(master.toString());
//        out.printf("Test diff objects%n");
//        out.println(test.toString());

        List<Difference> diffs = diffEngine.diff(master, test);
        if ( showItemizedDifferences ) {
            out.printf("Itemized results%n");
            for ( Difference diff : diffs )
                out.printf("DIFF: %s%n", diff.toString());
        }

        DiffEngine.SummaryReportParams params = new DiffEngine.SummaryReportParams(out, MAX_DIFFS, MAX_COUNT1_DIFFS, minCountForDiff);
        params.setDescending(false);
        diffEngine.reportSummarizedDifferences(diffs, params);
    }
}
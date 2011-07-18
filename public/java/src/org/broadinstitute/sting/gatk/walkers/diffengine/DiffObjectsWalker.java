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
 * Compares two record-oriented files, itemizing specific difference between equivalent
 * records in the two files.  Reports both itemized and summarized differences.
 * @author Mark DePristo
 * @version 0.1
 */
@Requires(value={})
public class DiffObjectsWalker extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(fullName="maxObjectsToRead", shortName="motr", doc="Max. number of objects to read from the files.  -1 [default] means unlimited", required=false)
    int MAX_OBJECTS_TO_READ = -1;

    @Argument(fullName="maxDiffs", shortName="M", doc="Max. number of diffs to process", required=false)
    int MAX_DIFFS = 0;

    @Argument(fullName="maxCount1Diffs", shortName="M1", doc="Max. number of diffs occuring exactly once in the file to process", required=false)
    int MAX_COUNT1_DIFFS = 0;

    @Argument(fullName="minCountForDiff", shortName="MCFD", doc="Min number of observations for a records to display", required=false)
    int minCountForDiff = 1;

    @Argument(fullName="showItemizedDifferences", shortName="SID", doc="Should we enumerate all differences between the files?", required=false)
    boolean showItemizedDifferences = false;

    @Argument(fullName="master", shortName="m", doc="Master file: expected results", required=true)
    File masterFile;

    @Argument(fullName="test", shortName="t", doc="Test file: new results to compare to the master file", required=true)
    File testFile;

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
        diffEngine.reportSummarizedDifferences(diffs, params);
    }
}
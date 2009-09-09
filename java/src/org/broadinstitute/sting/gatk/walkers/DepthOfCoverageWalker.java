/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import net.sf.samtools.SAMRecord;

/**
 * Display the depth of coverage at a given locus.
 */
public class DepthOfCoverageWalker extends LocusWalker<Integer, Pair<Long, Long>> {
    enum printType {
        NONE,
        COMPACT,
        DETAILED
    }

    @Argument(fullName="printStyle", shortName = "s", doc="Printing style: NONE, COMPACT, or DETAILED", required=false)
    public printType printStyle = printType.COMPACT;

    @Argument(fullName="excludeDeletions", shortName = "ed", doc="If true, we will exclude reads with deletions at a locus in coverage",required=false)
    public boolean excludeDeletionsInCoverage = false;

    @Argument(fullName="minMAPQ", shortName = "minMAPQ", doc="If provided, we will exclude reads with MAPQ < this value at a locus in coverage",required=false)
    public int excludeMAPQBelowThis = -1;

    public boolean includeReadsWithDeletionAtLoci() { return ! excludeDeletionsInCoverage; }

    public void initialize() {
        switch ( printStyle ) {
            case COMPACT:
                out.printf("locus depth%n");
                break;
            case DETAILED:
                out.printf("locus nCleanReads nDeletionReads nLowMAPQReads%n");
                break;
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int nCleanReads = 0, nBadMAPQReads = 0, nDeletionReads = 0;

        for ( int i = 0; i < context.getReads().size(); i++ ) {
            SAMRecord read = context.getReads().get(i);
            int offset = context.getOffsets().get(i);

            if ( read.getMappingQuality() < excludeMAPQBelowThis ) nBadMAPQReads++;
            else if ( offset == -1 ) nDeletionReads++;
            else nCleanReads++;
        }

        int nTotalReads = nCleanReads + (excludeDeletionsInCoverage ? 0 : nDeletionReads);

        switch ( printStyle ) {
            case COMPACT:
                out.printf("%s %8d%n", context.getLocation(), nTotalReads);
                break;
            case DETAILED:
                out.printf("%s %8d %8d %8d %8d%n", context.getLocation(), nTotalReads, nCleanReads, nDeletionReads, nBadMAPQReads);
                break;
        }

        return nTotalReads;
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public Pair<Long, Long> reduceInit() { return new Pair<Long,Long>(0l,0l); }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair<Long,Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n", 
                ((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
    }
}
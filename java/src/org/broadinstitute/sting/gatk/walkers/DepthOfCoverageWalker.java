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
import org.broadinstitute.sting.utils.pileup.*;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.*;

/**
 * Display the depth of coverage at a given locus.
 */
@By(DataSource.REFERENCE)
public class DepthOfCoverageWalker extends LocusWalker<Integer, Pair<Long, Long>> {
    enum printType {
        NONE,
        COMPACT,
        DETAILED
    }

    @Argument(fullName="printStyle", shortName = "s", doc="Printing style: NONE, COMPACT, or DETAILED", required=false)
    protected printType printStyle = printType.COMPACT;

    @Argument(fullName="excludeDeletions", shortName = "ed", doc="If true, we will exclude reads with deletions at a locus in coverage",required=false)
    protected boolean excludeDeletionsInCoverage = false;

    @Argument(fullName="minMAPQ", shortName = "minMAPQ", doc="If provided, we will exclude reads with MAPQ < this value at a locus in coverage",required=false)
    protected int excludeMAPQBelowThis = -1;

    @Argument(fullName="byReadGroup", shortName="byRG", doc="List read depths for each read group")
    protected boolean byReadGroup = false;

    @Argument(fullName="bySample", shortName="bySample", doc="List read depths for each sample")
    protected boolean bySample = false;


    private TreeSet<String> readGroupNames = new TreeSet<String>();
    private TreeSet<String> sampleNames = new TreeSet<String>();



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

        if ( byReadGroup ) {
            List<SAMReadGroupRecord> readGroups = this.getToolkit().getSAMFileHeader().getReadGroups();
            for ( SAMReadGroupRecord record : readGroups )
                readGroupNames.add(record.getReadGroupId());
        }

        if ( bySample ) {
            List<SAMReadGroupRecord> readGroups = this.getToolkit().getSAMFileHeader().getReadGroups();
            for ( SAMReadGroupRecord record : readGroups ) {
                String sample = record.getSample();
                if ( sample != null )
                    sampleNames.add(sample);
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int nCleanReads = 0, nBadMAPQReads = 0, nDeletionReads = 0;

        HashMap<String, Integer> depthByReadGroup = new HashMap<String, Integer>();
        for ( String readGroupName : readGroupNames )
            depthByReadGroup.put(readGroupName, 0);
        HashMap<String, Integer> depthBySample = new HashMap<String, Integer>();
        for ( String sample : sampleNames )
            depthBySample.put(sample, 0);

        ReadBackedPileup pileup = context.getPileup();
        for ( PileupElement p : pileup ) {

            if ( p.getRead().getMappingQuality() < excludeMAPQBelowThis ) nBadMAPQReads++;
            else if ( p.isDeletion() ) nDeletionReads++;
            else nCleanReads++;

            SAMReadGroupRecord readGroup = p.getRead().getReadGroup();
            if ( readGroup == null )
                continue;

            if ( byReadGroup ) {
                String readGroupName = readGroup.getReadGroupId();
                int oldDepth = depthByReadGroup.get(readGroupName);
                depthByReadGroup.put(readGroupName, oldDepth + 1);
            }

            if ( bySample ) {
                String sample = readGroup.getSample();
                if ( sample != null ) {
                    int oldDepth = depthBySample.get(sample);
                    depthBySample.put(sample, oldDepth + 1);
                }
            }
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

        if ( byReadGroup ) {
            for ( String rg : readGroupNames ) {
                out.printf("  %s %8d%n", rg, depthByReadGroup.get(rg));
            }
        }

        if ( bySample ) {
            for ( String sample : sampleNames ) {
                out.printf("  %s %8d%n", sample, depthBySample.get(sample));
            }
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
                ((double)result.getFirst() / result.getSecond()), result.getFirst(), result.getSecond());
    }
}

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
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.*;

/**
 * Display the depth of coverage at a given locus.
 */
@By(DataSource.REFERENCE)
public class DepthOfCoverageWalker extends LocusWalker<DepthOfCoverageWalker.DoCInfo, DepthOfCoverageWalker.DoCInfo> {
    enum printType {
        NONE,
        COMPACT,
        DETAILED
    }

    @Argument(fullName="printStyle", shortName = "s", doc="Printing style: NONE, COMPACT, or DETAILED", required=false)
    protected printType printStyle = printType.COMPACT;

    @Argument(fullName="minMAPQ", shortName = "minMAPQ", doc="If provided, we will also list read counts with MAPQ >= this value at a locus in coverage",required=false)
    protected int excludeMAPQBelowThis = -1;

    @Argument(fullName="minDepth", shortName = "minDepth", doc="If provided, we will also list the percentage of loci with depth >= this value per interval",required=false)
    protected int minDepthForPercentage = -1;

    @Argument(fullName="byReadGroup", shortName="byRG", doc="List read depths for each read group")
    protected boolean byReadGroup = false;

    @Argument(fullName="bySample", shortName="bySample", doc="List read depths for each sample")
    protected boolean bySample = false;

    // keep track of the read group and sample names
    private TreeSet<String> readGroupNames = new TreeSet<String>();
    private TreeSet<String> sampleNames = new TreeSet<String>();

    public boolean includeReadsWithDeletionAtLoci() { return true; }

    public void initialize() {

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

        StringBuilder header = new StringBuilder("location\ttotal_coverage\tcoverage_without_deletions");
        if ( excludeMAPQBelowThis > 0 ) {
            header.append("\tcoverage_atleast_MQ");
            header.append(excludeMAPQBelowThis);
        }
        if ( byReadGroup ) {
            for ( String rg : readGroupNames ) {
                header.append("\tcoverage_for_");
                header.append(rg);
            }
        }
        if ( bySample ) {
            for ( String sample : sampleNames ) {
                header.append("\tcoverage_for_");
                header.append(sample);
            }
        }
        out.println(header.toString());
    }

    public DoCInfo map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        ReadBackedPileup pileup = context.getPileup();

        DoCInfo info = new DoCInfo();
        info.totalCoverage = pileup.size();

        int nBadMAPQReads = 0, nDeletionReads = 0;
        for ( PileupElement p : pileup ) {

            if ( excludeMAPQBelowThis > 0 && p.getRead().getMappingQuality() < excludeMAPQBelowThis )
                nBadMAPQReads++;
            else if ( p.isDeletion() )
                nDeletionReads++;

            SAMReadGroupRecord readGroup = p.getRead().getReadGroup();
            if ( readGroup == null )
                continue;

            if ( byReadGroup ) {
                String readGroupName = readGroup.getReadGroupId();
                int oldDepth = info.depthByReadGroup.get(readGroupName);
                info.depthByReadGroup.put(readGroupName, oldDepth + 1);
            }

            if ( bySample ) {
                String sample = readGroup.getSample();
                if ( sample != null ) {
                    int oldDepth = info.depthBySample.get(sample);
                    info.depthBySample.put(sample, oldDepth + 1);
                }
            }
        }

        info.numDeletions = nDeletionReads;
        if ( excludeMAPQBelowThis > 0 )
            info.numBadMQReads = nBadMAPQReads;

        printDoCInfo(context.getLocation(), info, false);

        return info;
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public DoCInfo reduceInit() { return new DoCInfo(); }

    public DoCInfo reduce(DoCInfo value, DoCInfo sum) {
        sum.totalCoverage += value.totalCoverage;
        sum.numDeletions += value.numDeletions;
        sum.numBadMQReads += value.numBadMQReads;
        if ( value.totalCoverage >= minDepthForPercentage ) {
            sum.minDepthCoveredLoci++;
        }
        if ( byReadGroup ) {
            for ( String rg : readGroupNames ) {
                int oldDepth = sum.depthByReadGroup.get(rg);
                sum.depthByReadGroup.put(rg, oldDepth + value.depthByReadGroup.get(rg));
            }
        }
        if ( bySample ) {
            for ( String sample : sampleNames ) {
                int oldDepth = sum.depthBySample.get(sample);
                sum.depthBySample.put(sample, oldDepth + value.depthBySample.get(sample));
            }
        }

        return sum;
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, DoCInfo>> results) {

        StringBuilder header = new StringBuilder("\nlocation\ttotal_coverage\taverage_coverage\tcoverage_without_deletions\taverage_coverage_without_deletions");
        if ( excludeMAPQBelowThis > 0 ) {
            header.append("\tcoverage_atleast_MQ");
            header.append(excludeMAPQBelowThis);
            header.append("\taverage_coverage_atleast_MQ");
            header.append(excludeMAPQBelowThis);
        }
        if ( minDepthForPercentage >= 0 ) {
            header.append("\tpercent_loci_covered_atleast_depth");
            header.append(minDepthForPercentage);            
        }
        if ( byReadGroup ) {
            for ( String rg : readGroupNames ) {
                header.append("\tcoverage_for_");
                header.append(rg);
            }
        }
        if ( bySample ) {
            for ( String sample : sampleNames ) {
                header.append("\tcoverage_for_");
                header.append(sample);
            }
        }
        out.println(header.toString());

        for ( Pair<GenomeLoc, DoCInfo> result : results )
            printDoCInfo(result.first, result.second, true);
    }

    private void printDoCInfo(GenomeLoc loc, DoCInfo info, boolean printAverageCoverage) {

        double totalBases = (double)(loc.getStop() - loc.getStart() + 1);

        StringBuilder sb = new StringBuilder();
        sb.append(loc);
        sb.append("\t");
        sb.append(info.totalCoverage);
        sb.append("\t");
        if ( printAverageCoverage ) {
            sb.append(String.format("%.2f", ((double)info.totalCoverage) / totalBases));
            sb.append("\t");
        }
        sb.append((info.totalCoverage - info.numDeletions));
        if ( printAverageCoverage ) {
            sb.append("\t");
            sb.append(String.format("%.2f", ((double)(info.totalCoverage - info.numDeletions)) / totalBases));
        }
        if ( excludeMAPQBelowThis > 0 ) {
            sb.append("\t");
            sb.append((info.totalCoverage - info.numBadMQReads));
            if ( printAverageCoverage ) {
                sb.append("\t");
                sb.append(String.format("%.2f", ((double)(info.totalCoverage - info.numBadMQReads)) / totalBases));
            }
        }
        if ( minDepthForPercentage >= 0 ) {
            sb.append("\t");
            sb.append(String.format("%.2f", ((double)info.minDepthCoveredLoci) / totalBases));            
        }
        if ( byReadGroup ) {
            for ( String rg : readGroupNames ) {
                sb.append("\t");
                sb.append(String.format("%8d", info.depthByReadGroup.get(rg)));
            }
        }

        if ( bySample ) {
            for ( String sample : sampleNames ) {
                sb.append("\t");
                sb.append(String.format("%8d", info.depthBySample.get(sample)));
            }
        }

        out.println(sb.toString());
    }

    public class DoCInfo {
        public int totalCoverage = 0;
        public int numDeletions = 0;
        public int numBadMQReads = 0;
        public int minDepthCoveredLoci = 0;

        public HashMap<String, Integer> depthByReadGroup = null;
        public HashMap<String, Integer> depthBySample = null;

        public DoCInfo() {
            if ( byReadGroup ) {
                depthByReadGroup = new HashMap<String, Integer>();
                for ( String readGroupName : readGroupNames )
                    depthByReadGroup.put(readGroupName, 0);
            }
            if ( bySample ) {
                depthBySample = new HashMap<String, Integer>();
                for ( String sample : sampleNames )
                    depthBySample.put(sample, 0);
            }
        }

    }
}

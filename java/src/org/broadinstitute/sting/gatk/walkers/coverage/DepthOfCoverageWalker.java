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

package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.gatk.walkers.*;
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

    @Argument(fullName="printBaseCounts", shortName = "bases", doc="Print individual base counts (A,C,G,T only)", required=false)
    protected boolean printBaseCounts = false;

    @Argument(fullName="minMAPQ", shortName = "minMAPQ", doc="If provided, we will also list read counts with MAPQ >= this value at a locus in coverage",required=false)
    protected int excludeMAPQBelowThis = -1;

    @Argument(fullName="minDepth", shortName = "minDepth", doc="If provided, we will also list the percentage of loci with depth >= this value per interval",required=false)
    protected int minDepthForPercentage = -1;

    @Argument(fullName="byReadGroup", shortName="byRG", doc="List read depths for each read group")
    protected boolean byReadGroup = false;

    @Argument(fullName="bySample", shortName="bySample", doc="List read depths for each sample")
    protected boolean bySample = false;

    @Argument(fullName="printHistogram", shortName="histogram", doc="Print a histogram of the coverage")
    protected boolean printHistogram = false;


    // keep track of the read group and sample names
    private TreeSet<String> readGroupNames = new TreeSet<String>();
    private TreeSet<String> sampleNames = new TreeSet<String>();

    // keep track of the histogram data
    private ExpandingArrayList<Integer> coverageHist = null;
    private int maxDepth = 0;
    private int totalLoci = 0;

    // we want to see reads with deletions
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    public void initialize() {

        // initialize histogram array
        if ( printHistogram ) {
            coverageHist = new ExpandingArrayList<Integer>();
        }

        // initialize read group names from BAM header
        if ( byReadGroup ) {
            List<SAMReadGroupRecord> readGroups = this.getToolkit().getSAMFileHeader().getReadGroups();
            for ( SAMReadGroupRecord record : readGroups )
                readGroupNames.add(record.getReadGroupId());
        }

        // initialize sample names from BAM header
        if ( bySample ) {
            List<SAMReadGroupRecord> readGroups = this.getToolkit().getSAMFileHeader().getReadGroups();
            for ( SAMReadGroupRecord record : readGroups ) {
                String sample = record.getSample();
                if ( sample != null )
                    sampleNames.add(sample);
            }
        }

        // build and print the per-locus header
        out.println("\nPER_LOCUS_COVERAGE_SECTION");
        StringBuilder header = new StringBuilder("location\ttotal_coverage\tcoverage_without_deletions");
        if ( excludeMAPQBelowThis > 0 ) {
            header.append("\tcoverage_atleast_MQ");
            header.append(excludeMAPQBelowThis);
        }
        if ( printBaseCounts ) {
            header.append("\tA_count\tC_count\tG_count\tT_count");
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

        // fill in and print all of the per-locus coverage data, then return it to reduce

        ReadBackedPileup pileup = context.getPileup();

        DoCInfo info = new DoCInfo();
        info.totalCoverage = pileup.size();

        int nBadMAPQReads = 0, nDeletionReads = 0;
        for ( PileupElement p : pileup ) {

            if ( excludeMAPQBelowThis > 0 && p.getRead().getMappingQuality() < excludeMAPQBelowThis )
                nBadMAPQReads++;
            else if ( p.isDeletion() )
                nDeletionReads++;

            if ( printBaseCounts ) {
                int baseIndex = BaseUtils.simpleBaseToBaseIndex(p.getBase());
                if ( baseIndex != -1 )
                    info.baseCounts[baseIndex]++;
            }

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

        // if we need to print the histogram, fill in the data
        if ( printHistogram )
            incCov(info.totalCoverage);

        printDoCInfo(context.getLocation(), info, false);

        return info;
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public DoCInfo reduceInit() { return new DoCInfo(); }

    public DoCInfo reduce(DoCInfo value, DoCInfo sum) {

        // combine all of the per-locus data for a given interval

        sum.totalCoverage += value.totalCoverage;
        sum.numDeletions += value.numDeletions;
        sum.numBadMQReads += value.numBadMQReads;
        if ( value.totalCoverage >= minDepthForPercentage ) {
            sum.minDepthCoveredLoci++;
        }
        if ( printBaseCounts ) {
            for (int baseIndex = 0; baseIndex < BaseUtils.BASES.length; baseIndex++ )
                sum.baseCounts[baseIndex] += value.baseCounts[baseIndex];
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

        // build and print the per-interval header
        out.println("\n\nPER_INTERVAL_COVERAGE_SECTION");

        StringBuilder header = new StringBuilder("location\ttotal_coverage\taverage_coverage\tcoverage_without_deletions\taverage_coverage_without_deletions");
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
        if ( printBaseCounts ) {
            header.append("\tA_count\tC_count\tG_count\tT_count");
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

        // print all of the individual per-interval coverage data
        for ( Pair<GenomeLoc, DoCInfo> result : results )
            printDoCInfo(result.first, result.second, true);

        // if we need to print the histogram, do so now
        if ( printHistogram )
            printHisto();
    }

    private void incCov(int depth) {
        int c = coverageHist.expandingGet(depth, 0);
        coverageHist.set(depth, c + 1);
        if ( depth > maxDepth )
            maxDepth = depth;
        totalLoci++;
    }

    private int getCov(int depth) {
        return coverageHist.get(depth);
    }

    private void printHisto() {

        // sanity check
        if ( totalLoci == 0 )
            return;

        // Code for calculting std devs adapted from Michael Melgar's python script

        // Find the maximum extent of 'good' data
        // First, find the mode
        long maxValue = getCov(1); // ignore doc=0
        int mode = 1;
        for (int i = 2; i <= maxDepth; i++) {
            if ( getCov(i) > maxValue ) {
                maxValue = getCov(i);
                mode = i;
            }
        }

        // now, procede to find end of good Gaussian fit
        long dist = (long)Math.pow(10, 9);
        while ( Math.abs(getCov(mode) - getCov(1)) < dist && mode < maxDepth )
            dist = Math.abs(getCov(mode++) - getCov(1));
        int maxGoodDepth = Math.min(mode + 1, maxDepth);

        // calculate the mean of the good region
        long totalGoodSites = 0, totalGoodDepth = 0;
        for (int i = 1; i <= maxGoodDepth; i++) { // ignore doc=0
            totalGoodSites += getCov(i);
            totalGoodDepth += i * getCov(i);
        }
        double meanGoodDepth = (double)totalGoodDepth / (double)totalGoodSites;

        // calculate the variance and standard deviation of the good region
        double var = 0.0;
        for (int i = 1; i <= maxGoodDepth; i++) {  // ignore doc=0
            var += getCov(i) * Math.pow(meanGoodDepth - (double)i, 2);
        }
        double stdev = Math.sqrt(var / (double)totalGoodSites);

        // print
        out.println("\n\nHISTOGRAM_SECTION");
        out.printf("# sites within Gaussian fit  : mean:%f num_sites:%d std_dev:%f%n", meanGoodDepth, totalGoodSites, stdev);

        for (int i = 1; i <= 5; i++)
            out.printf("# Gaussian mean + %d Std Dev  : %f%n", i, (meanGoodDepth + i*stdev));

		out.println("\ndepth count freq(percent)");
		for (int i = 0; i <= maxDepth; i++)
			out.printf("%d %d %f\n", i, getCov(i), (100.0*getCov(i)) / (double)totalLoci);
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

        if ( printBaseCounts ) {
            for (int baseIndex = 0; baseIndex < BaseUtils.BASES.length; baseIndex++ ) {
                sb.append("\t");
                sb.append(String.format("%8d", info.baseCounts[baseIndex]));
            }
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

        public int[] baseCounts = null;

        public HashMap<String, Integer> depthByReadGroup = null;
        public HashMap<String, Integer> depthBySample = null;

        public DoCInfo() {
            if ( printBaseCounts ) {
                baseCounts = new int[4];
            }
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

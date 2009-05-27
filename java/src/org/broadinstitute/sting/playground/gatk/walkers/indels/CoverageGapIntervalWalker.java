package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

@WalkerName("CoverageGapIntervals")
public class CoverageGapIntervalWalker extends LocusWalker<Pair<GenomeLoc, Integer>, GenomeLoc> {

    private final int minReadsAtInterval = 10;

    public void initialize() {}

    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        int goodReads = 0;
        List<SAMRecord> reads = context.getReads();
        for (int i = 0; i < reads.size(); i++ ) {
            if ( reads.get(i).getMappingQuality() > 0 )
                goodReads++;
        }
        return goodReads >= minReadsAtInterval;
    }

     public Pair<GenomeLoc, Integer> map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        // find the probability that this locus has a statistically significant gap in coverage
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        int totalXi = 0;
        for (int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            if ( read.getMappingQuality() == 0 )
                continue;
            int halfLength = read.getReadString().length() >> 1;
            int distanceFromMiddle = Math.abs(offsets.get(i) - halfLength);
            int quarterLength = halfLength >> 1;

            // Xi is < 0 if you are closer to the middle than the quartile
            //  and is > 0 if further to the middle than quartile
            // We expect the total sum of Xi over an interval to be ~0
            int Xi = distanceFromMiddle - quarterLength;
            totalXi += Xi;
        }

        return new Pair<GenomeLoc, Integer>(context.getLocation(), totalXi);
    }

    public void onTraversalDone() {}

    public GenomeLoc reduceInit() {
        return null;
    }

    public GenomeLoc reduce(Pair<GenomeLoc, Integer> value, GenomeLoc sum) {
        if ( value.second > 1000 ) {
            if ( sum != null )
                sum.setStop(value.first.getStop());
            else
                sum = value.first;
        } else if ( sum != null ) {
            out.println(sum);
            sum = null;
        }
        return sum;
    }
}
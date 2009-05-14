
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

@WalkerName("EntropyIntervals")
public class EntropyIntervalWalker extends LocusWalker<Pair<GenomeLoc, Double>, Pair<LinkedList<Double>, GenomeLoc>> {
    @Argument(fullName="windowSize", shortName="window", doc="window size for calculating entropy", required=false)
    public int windowSize = 10;

    public void initialize() {
        if ( windowSize < 1)
            throw new RuntimeException("Window Size must be a positive integer");
    }

    public Pair<GenomeLoc, Double> map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        // return the entropy of this locus
        int[] baseCounts = new int[4];
        for (int i=0; i < 4; i++)
            baseCounts[i] = 0;

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        int goodBases = 0;
        double errorRate = 0.0;
        for (int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            int base = BaseUtils.simpleBaseToBaseIndex((char)read.getReadBases()[offset]);
            if ( base != -1 ) {
                errorRate += Math.pow(10.0, (double)read.getBaseQualities()[offset] / -10.0);
                baseCounts[base]++;
                goodBases++;
            }
        }

        double expectedEntropy = (errorRate * Math.log(errorRate)) + ((1-errorRate) * Math.log(1-errorRate));
        double observedEntropy = 0.0;
        if ( goodBases > 0 ) {
            for (int i=0; i < 4; i++) {
                double Pjk = (double)baseCounts[i] / (double)goodBases;
                if ( Pjk > 0 )
                    observedEntropy += Pjk * Math.log(Pjk);
            }
            if ( observedEntropy != 0 )
                observedEntropy *= -1;
        }

        double locusEntropy = (observedEntropy > expectedEntropy ? (observedEntropy-expectedEntropy) : 0.0);
        return new Pair<GenomeLoc, Double>(context.getLocation(), locusEntropy);
    }

    public void onTraversalDone() {}

    public Pair<LinkedList<Double>, GenomeLoc> reduceInit() {
        return new Pair<LinkedList<Double>, GenomeLoc>(new LinkedList<Double>(), null);
    }

    public Pair<LinkedList<Double>, GenomeLoc> reduce(Pair<GenomeLoc, Double> value, Pair<LinkedList<Double>, GenomeLoc> sum) {
        sum.first.addLast(value.second);
        if ( sum.first.size() <= windowSize )
            return sum;

        sum.first.remove();
        double avgEntropy = 0.0;
        for (int i = 0; i < windowSize; i++)
            avgEntropy += sum.first.get(i);
        avgEntropy /= windowSize;

        if ( avgEntropy > 0.001 ) {
            //out.println(avgEntropy);

            // if there is no interval to the left, then this is the first one
            if ( sum.second == null ) {
                sum.second = value.first;
            }
            // if the intervals don't overlap, print out the leftmost one and start a new one
            else if ( !sum.second.contiguousP(value.first) ) {
                out.println(sum.second);
                sum.second = value.first;
            }
            // otherwise, merge them
            else {
                sum.second = sum.second.merge(value.first);
            }
        }

        return sum;
    }
}
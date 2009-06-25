
package org.broadinstitute.sting.playground.gatk.walkers.indels;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

@WalkerName("MismatchIntervals")
public class MismatchIntervalWalker extends LocusWalker<Pair<GenomeLoc, Boolean>, Pair<LinkedList<Boolean>, GenomeLoc>> {
    @Argument(fullName="windowSize", shortName="window", doc="window size for calculating entropy", required=false)
    public int windowSize = 10;
    @Argument(fullName="mismatchFraction", shortName="mismatch", doc="fraction of mismatching base qualities threshold", required=false)
    public double mismatchThreshold = 0.15;
    @Argument(fullName="allow454Reads", shortName="454", doc="process 454 reads", required=false)
    public boolean allow454 = false;

    private final int minReadsAtInterval = 4;

    public void initialize() {
        if ( windowSize < 1)
            throw new RuntimeException("Window Size must be a positive integer");
    }

    public Pair<GenomeLoc, Boolean> map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        char upperRef = Character.toUpperCase(ref);
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        int goodReads = 0, mismatchQualities = 0, totalQualities = 0;
        for (int i = 0; i < reads.size(); i++) {
            SAMRecord read = reads.get(i);
            if ( read.getMappingQuality() == 0 ||
                 read.getAlignmentBlocks().size() > 1 ||
		 (!allow454 && Utils.is454Read(read, getToolkit().getEngine().getSAMHeader())) )
                 continue;

            goodReads++;
            int offset = offsets.get(i);
            int quality = (int)read.getBaseQualityString().charAt(offset) - 33;
            totalQualities += quality;

            char base = Character.toUpperCase((char)read.getReadBases()[offset]);
            if ( base != upperRef )
                mismatchQualities += quality;
        }

        boolean flag = false;
        if ( goodReads >= minReadsAtInterval && (double)mismatchQualities / (double)totalQualities > mismatchThreshold )
            flag = true;
        return new Pair<GenomeLoc, Boolean>(context.getLocation(), flag);
    }

    public void onTraversalDone(Pair<LinkedList<Boolean>, GenomeLoc> sum) {
        if (sum.second != null)
            out.println(sum.second);
    }

    public Pair<LinkedList<Boolean>, GenomeLoc> reduceInit() {
        return new Pair<LinkedList<Boolean>, GenomeLoc>(new LinkedList<Boolean>(), null);
    }

    public Pair<LinkedList<Boolean>, GenomeLoc> reduce(Pair<GenomeLoc, Boolean> value, Pair<LinkedList<Boolean>, GenomeLoc> sum) {
        // if we hit a new contig, clear the list
        if ( sum.second != null && sum.second.getContigIndex() != value.first.getContigIndex() ) {
            sum.first.clear();
            sum.second = null;
        }
        
        sum.first.addLast(value.second);
        if ( sum.first.size() <= windowSize )
            return sum;

        sum.first.remove();
        if ( !value.second )
            return sum;

        int mismatches = 0;
        int firstMismatch = -1;
        for (int i = 0; i < windowSize; i++) {
            if ( sum.first.get(i) ) {
                mismatches++;
                if ( firstMismatch == -1 )
                    firstMismatch = i;
            }
        }

        if ( mismatches > 1 ) {
            // if there is no interval to the left, then this is the first one
            if ( sum.second == null ) {
                sum.second = value.first;
                sum.second.setStart(sum.second.getStart() - windowSize + firstMismatch + 1);
            }
            // if the intervals don't overlap, print out the leftmost one and start a new one
            else if ( value.first.getStop() - sum.second.getStop() > windowSize ) {
                out.println(sum.second);
                sum.second = value.first;
                sum.second.setStart(sum.second.getStart() - windowSize + firstMismatch + 1);
            }
            // otherwise, merge them
            else {
                sum.second.setStop(value.first.getStop());
            }
        }

        return sum;
    }
}

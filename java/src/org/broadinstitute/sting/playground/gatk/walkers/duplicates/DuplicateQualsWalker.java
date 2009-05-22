package org.broadinstitute.sting.playground.gatk.duplicates.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.duplicates.DuplicateComp;
import org.broadinstitute.sting.utils.duplicates.DupUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.List;
import java.util.ArrayList;
import java.io.PrintStream;

import net.sf.samtools.SAMRecord;

class MismatchCounter {
    long nObs = 0;
    long nMismatches = 0;

    public void inc(long incNObs, long incNMismatches) {
        nObs += incNObs;
        nMismatches += incNMismatches;
    }

    public void inc(boolean mismatchP) {
        inc(1, mismatchP ? 1 : 0);
    }


    public double mismatchRate() {
        return (double)nMismatches / nObs;
    }

    public byte empiricalQualScore() {
        return QualityUtils.probToQual(1 - mismatchRate(), 0);
    }

    public String headerString() {
        return "mismatchRate\tempiricalQ\tnObs\tnMismatches";
    }

    public String toString() {
        return String.format("%.10f\t%d\t%d\t%6d", mismatchRate(), empiricalQualScore(), nObs, nMismatches);
    }
}

class QualityTracker {
    final private int MAX_QUAL_SCORE = 100;
    MismatchCounter[][] mismatchesByQ = new MismatchCounter[MAX_QUAL_SCORE][MAX_QUAL_SCORE];

    public QualityTracker() {
        for ( int i = 0; i < MAX_QUAL_SCORE; i++ ) {
            for ( int j = 0; j < MAX_QUAL_SCORE; j++ ) {
                mismatchesByQ[i][j] = new MismatchCounter();
            }
        }
    }

    public void inc(int b1Qi, int b2Qi, boolean mismatchP, boolean orderDependent) {
        int b1Q = orderDependent ? b1Qi : Math.max(b1Qi, b2Qi);
        int b2Q = orderDependent ? b2Qi : Math.min(b1Qi, b2Qi);

        if ( b1Q > MAX_QUAL_SCORE ) throw new RuntimeException("Unexpectedly large base quality " + b1Q);
        if ( b2Q > MAX_QUAL_SCORE ) throw new RuntimeException("Unexpectedly large base quality " + b2Q);

        mismatchesByQ[b1Q][b2Q].inc(mismatchP);
    }

    public void inc(DuplicateComp dc, boolean orderDependent) {
        inc(dc.getQLarger(), dc.getQSmaller(), dc.isMismatchP(), orderDependent);
    }

    public int probMismatchQ1Q2(int q1, int q2) {
        double e1 = 1 - QualityUtils.qualToProb(q1);
        double e2 = 1 - QualityUtils.qualToProb(q2);
        double eMM = e1 * (1 - e2) + (1 - e1) * e2 - 1/3 * e1 * e2;
        return QualityUtils.probToQual(1 - eMM, 0.0);
    }
    
    public void printToStream(PrintStream out, boolean filterUnobserved) {
        out.printf("Q1\tQ2\tQmin\t%s%n", mismatchesByQ[0][0].headerString());
        for ( int i = 0; i < MAX_QUAL_SCORE; i++ ) {
            for ( int j = 0; j < MAX_QUAL_SCORE; j++ ) {
                MismatchCounter mc = mismatchesByQ[i][j];
                //System.out.printf("MC = %s%n", mc);
                if ( filterUnobserved && mc.nObs == 0 )
                    continue;
                out.printf("%d\t%d\t%d\t%s\t%n", i, j, probMismatchQ1Q2(i,j), mc.toString());
            }
        }
    }
}

public class DuplicateQualsWalker extends DuplicateWalker<List<DuplicateComp>, QualityTracker> {
    @Argument(fullName="filterUnobservedQuals", required=false, doc="Show only quality bins with at least one observation in the data")
    public boolean FILTER_UNOBSERVED_QUALS = false;

    @Argument(fullName="maxPairwiseCompsPerDupSet", required=false, doc="Maximumize number of pairwise comparisons to perform among duplicate read sets")
    public int MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET = 100;

    @Argument(fullName="combinedQuals", required=false, doc="Combine and assess pairwise base qualities")
    public boolean COMBINE_QUALS = false;

    @Argument(fullName="combineAllDups", required=false, doc="Combine and assess pairwise base qualities")
    public boolean COMBINE_ALL_DUPS = false;

    @Argument(fullName="orderDependent", required=false, doc="")
    public boolean orderDependent = false;

    @Argument(fullName="compareToUniqueReads", required=false, doc="If true, then we will compare only to unique (i.e., non-duplicated molecules) at the same duplicate site")
    public boolean compareToUniqueReads = false;

    @Argument(fullName="comparePairToSingleton", required=false, doc="If true, then we will compare a combined dup to a random other read in the duplicate set, not a combined pair itself")
    public boolean comparePairToSingleton = false;

    final boolean DEBUG = false;
    final private boolean ACTUALLY_DO_WORK = true;

    public void onTraversalDone(QualityTracker result) {
        result.printToStream(out, FILTER_UNOBSERVED_QUALS);
    }

    public QualityTracker reduceInit() {
        return new QualityTracker();
    }

    public QualityTracker reduce(List<DuplicateComp> dupComps, QualityTracker tracker) {
        for ( DuplicateComp dc : dupComps ) {
            tracker.inc(dc, orderDependent);
        }
        
        return tracker;
    }

    // Print out data for regression
    public List<DuplicateComp> map(GenomeLoc loc, byte[] refBases, LocusContext context,
                                   List<SAMRecord> uniqueReads,
                                   List<SAMRecord> duplicateReads) {
        //logger.info(String.format("%s has %d duplicates and %d non-duplicates", loc, duplicateReads.size(), uniqueReads.size()));
        List<DuplicateComp> pairwiseComps = new ArrayList<DuplicateComp>();
        
        if ( ! ACTUALLY_DO_WORK )
            return pairwiseComps;

        if ( COMBINE_QUALS ) {
            Pair<SAMRecord, SAMRecord> combinedReads = DupUtils.combinedReadPair( duplicateReads );
            if ( combinedReads != null ) {
                SAMRecord combined1 = combinedReads.first;
                SAMRecord combined2 = combinedReads.second;

                if ( comparePairToSingleton )
                    pairwiseComps = addPairwiseMatches( pairwiseComps, combined1, duplicateReads.get(2), uniqueReads );
                else
                    pairwiseComps = addPairwiseMatches( pairwiseComps, combined1, combined2, uniqueReads );
            }
        } else {
            int nComparisons = 0;
            for ( SAMRecord read1 : duplicateReads ) {
                for ( SAMRecord read2 : duplicateReads ) {
                    if ( read1.hashCode() < read2.hashCode() && DupUtils.usableDuplicate(read1, read2) ) {
                        // the hashcode insures we don't do A vs. B and B vs. A
                        //System.out.printf("Comparing %s against %s%n", read1, read2);
                        nComparisons++;
                        pairwiseComps = addPairwiseMatches( pairwiseComps, read1, read2, uniqueReads );
                        if ( nComparisons > MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET )
                            break;
                    }
                }
            }
        }

        return pairwiseComps;
    }
    
    private List<DuplicateComp> addPairwiseMatches(List<DuplicateComp> comps,
                                                   SAMRecord read1, SAMRecord read2,
                                                   List<SAMRecord> uniqueReads ) {
        if ( compareToUniqueReads ) {
            // we want to compare to a read in the unique read set
            if ( uniqueReads.size() > 0 ) {                 // there's actually something to compare to
                SAMRecord uniqueRead = uniqueReads.get(0);  // might as well get the first one
                return pairwiseMatches(comps, read1, uniqueRead);
            } else {
                return comps;
            }
        } else {
            // default, just do read1 vs. read2
            return pairwiseMatches(comps, read1, read2);
        }
    }

    /**
     * Calculates the pairwise mismatches between reads read1 and read2 and adds the result to the comps list.
     * Doesn't contain any logic deciding what to compare, just does read1 and read2
     * 
     * @param comps
     * @param read1
     * @param read2
     * @return
     */
    private List<DuplicateComp> pairwiseMatches(List<DuplicateComp> comps, SAMRecord read1, SAMRecord read2 ) {
        byte[] read1Bases = read1.getReadBases();
        byte[] read1Quals = read1.getBaseQualities();
        byte[] read2Bases = read2.getReadBases();
        byte[] read2Quals = read2.getBaseQualities();

        for ( int i = 0; i < read1Bases.length; i++) {
            byte qual1 = read1Quals[i];
            byte qual2 = read2Quals[i];
            boolean mismatchP = ! BaseUtils.basesAreEqual(read1Bases[i], read2Bases[i]);
            DuplicateComp dc = new DuplicateComp(qual1, qual2, mismatchP);
            comps.add(dc);
        }

        return comps;
    }
}
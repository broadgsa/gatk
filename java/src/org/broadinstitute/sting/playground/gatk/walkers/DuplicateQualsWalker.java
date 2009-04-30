package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.List;
import java.util.ArrayList;
import java.io.PrintStream;
import java.io.FileNotFoundException;

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
    final private int MAX_QUAL_SCORE = 50;
    MismatchCounter[][] mismatchesByQ = new MismatchCounter[MAX_QUAL_SCORE][MAX_QUAL_SCORE];

    public QualityTracker() {
        for ( int i = 0; i < MAX_QUAL_SCORE; i++ ) {
            for ( int j = 0; j < MAX_QUAL_SCORE; j++ ) {
                mismatchesByQ[i][j] = new MismatchCounter();
            }
        }
    }

    public void inc(int b1Q, int b2Q, boolean mismatchP) {
        if ( b1Q > MAX_QUAL_SCORE ) throw new RuntimeException("Unexpectedly large base quality " + b1Q);
        if ( b2Q > MAX_QUAL_SCORE ) throw new RuntimeException("Unexpectedly large base quality " + b2Q);
        mismatchesByQ[b1Q][b2Q].inc(mismatchP);
    }

    public void inc(DuplicateComp dc) {
        inc(dc.qLarger, dc.qSmaller, dc.mismatchP);
    }

    public void printToStream(PrintStream out, boolean filterUnobserved) {
        out.printf("Q1\tQ2\t%s%n", mismatchesByQ[0][0].headerString());
        for ( int i = 0; i < MAX_QUAL_SCORE; i++ ) {
            for ( int j = 0; j < MAX_QUAL_SCORE; j++ ) {
                MismatchCounter mc = mismatchesByQ[i][j];
                //System.out.printf("MC = %s%n", mc);
                if ( filterUnobserved && mc.nObs == 0 )
                    continue;
                out.printf("%d\t%d\t%s\t%n", i, j, mc.toString());
            }
        }
    }
}

class DuplicateComp {
    int qLarger;
    int qSmaller;
    boolean mismatchP;

    public DuplicateComp( int b1Q, int b2Q, boolean mismatchP ) {
        qLarger = Math.max(b1Q, b2Q);
        qSmaller = Math.min(b1Q, b2Q);
        this.mismatchP = mismatchP;
    }

    public String toString() {
        return String.format("%d %d %b", qLarger, qSmaller, mismatchP);
    }
}

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class DuplicateQualsWalker extends DuplicateWalker<List<DuplicateComp>, QualityTracker> {
    @Argument(fullName="filterUnobservedQuals", shortName="filterUnobservedQuals", required=false, defaultValue="false", doc="Show only quality bins with at least one observation in the data")
    public boolean FILTER_UNOBSERVED_QUALS;

    @Argument(fullName="maxPairwiseCompsPerDupSet", shortName="maxPairwiseCompsPerDupSet", required=false, defaultValue="100", doc="Maximumize number of pairwise comparisons to perform among duplicate read sets")
    public int MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET;

    final private boolean ACTUALLY_DO_WORK = true;

    public void onTraversalDone(QualityTracker result) {
        result.printToStream(out, FILTER_UNOBSERVED_QUALS);
    }

    public QualityTracker reduceInit() {
        return new QualityTracker();
    }

    public QualityTracker reduce(List<DuplicateComp> dupComps, QualityTracker tracker) {
        for ( DuplicateComp dc : dupComps ) {
            tracker.inc(dc);
        }
        
        return tracker;
    }

    // Print out data for regression
    public List<DuplicateComp> map(GenomeLoc loc, byte[] refBases, LocusContext context, List<SAMRecord> duplicateReads) {
        //out.printf("%s has %d duplicates%n", loc, duplicateReads.size());
        List<DuplicateComp> comps = new ArrayList<DuplicateComp>();

        if ( ! ACTUALLY_DO_WORK )
            return comps;

        int nComparisons = 0;
        for ( SAMRecord read1 : duplicateReads ) {
            byte[] read1Bases = read1.getReadBases();
            byte[] read1Quals = read1.getBaseQualities();

            for ( SAMRecord read2 : duplicateReads ) {
                if ( read1 != read2 && read1.getReadLength() == read2.getReadLength()) {
                    byte[] read2Bases = read2.getReadBases();
                    byte[] read2Quals = read2.getBaseQualities();
                    nComparisons++;

                    for ( int i = 0; i < read1Bases.length; i++) {
                        byte read1Q = read1Quals[i];
                        byte read2Q = read2Quals[i];
                        boolean mismatchP = read1Bases[i] != read2Bases[i];
                        DuplicateComp dc = new DuplicateComp(read1Q, read2Q, mismatchP);

                        //logger.debug(String.format("dc: %s", dc));
                        comps.add(dc);
                    }

                    if ( nComparisons > MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET )
                        break;
                }
            }
        }

        return comps;
    }
}
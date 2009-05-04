package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;
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
    final private int MAX_QUAL_SCORE = 100;
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

    @Argument(fullName="combinedQuals", shortName="combinedQuals", required=false, defaultValue="false", doc="Combine and assess pairwise base qualities")
    public boolean COMBINE_QUALS;

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
            tracker.inc(dc);
        }
        
        return tracker;
    }

    // Print out data for regression
    public List<DuplicateComp> map(GenomeLoc loc, byte[] refBases, LocusContext context, List<SAMRecord> duplicateReads) {
        //out.printf("%s has %d duplicates%n", loc, duplicateReads.size());
        List<DuplicateComp> pairwiseComps = new ArrayList<DuplicateComp>();
        
        if ( ! ACTUALLY_DO_WORK )
            return pairwiseComps;

        if ( COMBINE_QUALS ) {
            Pair<SAMRecord, SAMRecord> combinedReads = combinedReadPair( duplicateReads );
            if ( combinedReads != null ) {
                SAMRecord combined1 = combinedReads.first;
                SAMRecord combined2 = combinedReads.second;
                addPairwiseMatches( pairwiseComps, combined1, combined2 );
            }
        }
        else {
            int nComparisons = 0;
            for ( SAMRecord read1 : duplicateReads ) {
                for ( SAMRecord read2 : duplicateReads ) {
                    if ( usableDuplicate(read1, read2) ) {
                        nComparisons++;
                        addPairwiseMatches( pairwiseComps, read1, read2 );
                        if ( nComparisons > MAX_PAIRSIZE_COMPS_PER_DUPLICATE_SET )
                            break;
                    }
                }
            }
        }

        return pairwiseComps;
    }

    private boolean usableDuplicate( SAMRecord read1, SAMRecord read2 ) {
            return read1 != read2 && read1.getReadLength() == read2.getReadLength();
    }
    
    private List<DuplicateComp> addPairwiseMatches(List<DuplicateComp> comps,
                                                 SAMRecord read1, SAMRecord read2 ) {
        byte[] read1Bases = read1.getReadBases();
        byte[] read1Quals = read1.getBaseQualities();
        byte[] read2Bases = read2.getReadBases();
        byte[] read2Quals = read2.getBaseQualities();

        for ( int i = 0; i < read1Bases.length; i++) {
            byte qual1 = read1Quals[i];
            byte qual2 = read2Quals[i];
            boolean mismatchP = read1Bases[i] != read2Bases[i];
            DuplicateComp dc = new DuplicateComp(qual1, qual2, mismatchP);
            comps.add(dc);
        }

        return comps;
    }

    private Pair<SAMRecord, SAMRecord> combinedReadPair( List<SAMRecord> duplicateReads ) {
        if ( duplicateReads.size() < 4 )
            return null;

        SAMRecord c1 = combineDuplicates(duplicateReads.get(0),duplicateReads.get(1));
        SAMRecord c2 = combineDuplicates(duplicateReads.get(2),duplicateReads.get(3));
        return new Pair<SAMRecord, SAMRecord>(c1, c2);
    }

    private SAMRecord sample3rdRead( List<SAMRecord> duplicateReads, SAMRecord read1, SAMRecord read2 ) {
         if ( duplicateReads.size() <= 2 ) {
             // no third unique read is available
             return null;
         } else {
             for ( SAMRecord read3 : duplicateReads ) {
                 if ( usableDuplicate(read1, read3) && usableDuplicate(read2, read3) )
                     return read3;
             }
             
             return null;
         }
    }

    private SAMRecord tmpCopyRead(SAMRecord read) {
        SAMRecord copy = new SAMRecord(read.getHeader());
        copy.setReadName(read.getReadName());
        //copy.setReadString(final String value) {
        copy.setReadBases(read.getReadBases());
        copy.setBaseQualities(read.getBaseQualities());
        copy.setReferenceName(read.getReferenceName());
        copy.setReferenceIndex(read.getReferenceIndex());
        copy.setMateReferenceName(read.getMateReferenceName());
        copy.setMateReferenceIndex(read.getMateReferenceIndex());
        copy.setAlignmentStart(read.getAlignmentStart());
                //copy.setAlignmentEnd(read.getAlignmentEnd());
        copy.setMateAlignmentStart(read.getMateAlignmentStart());
        copy.setInferredInsertSize(read.getInferredInsertSize());
        copy.setMappingQuality(read.getMappingQuality());
        copy.setCigar(read.getCigar());
        copy.setFlags(copy.getFlags());
        
        return copy;
    }

    private SAMRecord combineDuplicates(SAMRecord read1, SAMRecord read2) {
        byte[] read1Bases = read1.getReadBases();
        byte[] read1Quals = read1.getBaseQualities();
        byte[] read2Bases = read2.getReadBases();
        byte[] read2Quals = read2.getBaseQualities();

        byte[] bases = new byte[read1Bases.length];
        byte[] quals = new byte[read1Bases.length];

        SAMRecord c = tmpCopyRead(read1);
        for ( int i = 0; i < read1Bases.length; i++) {
            byte base1 = read1Bases[i];
            byte base2 = read2Bases[i];
            byte qual1 = read1Quals[i];
            byte qual2 = read2Quals[i];
            final double p1 = QualityUtils.qualToProb(qual1);
            final double p2 = QualityUtils.qualToProb(qual2);

            double pc;
            byte basec;

            if ( base1 == base2 ) {
                // agreement
                basec = base1;
                pc = 1 - (1 - p1) * (1 - p2);
            } else {
                // disagreement
                basec = p1 > p2 ? base1 : base2;
                pc    = p1 > p2 ? p1 : p2;
                //pc = 0;
            }

            bases[i] = basec;
            quals[i] = QualityUtils.probToQual(pc, 0.0);

            if ( DEBUG )
                logger.debug(String.format("Combining %s (Q%2d) with %s (Q%2d) -> %s (Q%2d)%s%n",
                    (char)base1, qual1, (char)base2, qual2, (char)bases[i], quals[i],
                    base1 == base2 ? "" : " [MISMATCH]"));
        }
        c.setReadBases(bases);
        c.setBaseQualities(quals);

        return c;
    }
}
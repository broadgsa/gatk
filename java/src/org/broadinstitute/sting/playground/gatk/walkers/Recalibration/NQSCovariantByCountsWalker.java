package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.List;
import java.util.LinkedList;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 29, 2009
 * Time: 10:54:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class NQSCovariantByCountsWalker extends LocusWalker< LocalMapType, int[][][][] > {

    final static int SMALL_DEVIATION = 5;
    final static int LARGE_DEVIATION = 8;
    final static int WIN_SIZE_SIDE = 4;
    final static int MATCH_OFFSET = 1;
    final static int MM_OFFSET = 0;
    final static int MAX_Q_SCORE = 2 + QualityUtils.MAX_REASONABLE_Q_SCORE;
    final static String DATA_FORMAT = "%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d%n";
    final static String HEADER_FORMAT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n";

    int NGHD_SIZE;

    public void initialize() {
        NGHD_SIZE = 2*WIN_SIZE_SIDE+1;
    }

    public int[][][][] reduceInit() {
        int[][][][] q =  new int[2*WIN_SIZE_SIDE][2*WIN_SIZE_SIDE][MAX_Q_SCORE][2];

        for ( int i = 0; i < 2*WIN_SIZE_SIDE; i ++ ) {
            for ( int j = 0; j < 2*WIN_SIZE_SIDE; j ++ ) {
                for ( int k = 0; k < MAX_Q_SCORE; k ++ ) {
                    for ( int h = 0; h < 2; h ++ ) {
                        q[i][j][k][h] = 0;
                    }
                }
            }
        }

        return q;
    }

    public LocalMapType map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        return new LocalMapType(context, ref, tracker);
    }

    public int [][][][] reduce( LocalMapType map, int [][][][] quals ) {
        for ( int i = 0; i < map.numReads(); i ++ ) {
            SAMRecord read = map.context.getReads().get(i);
            int offset = map.context.getOffsets().get(i);
            Pair<Integer,Integer> deviantCounts = countNumDeviantBases(read,offset);

            if ( isMismatch(read,offset,map.ref) ) {
                quals[deviantCounts.first][deviantCounts.second][(int) read.getBaseQualities()[offset]][MM_OFFSET]++;
            } else {
                quals[deviantCounts.first][deviantCounts.second][(int) read.getBaseQualities()[offset]][MATCH_OFFSET]++;
            }
        }

        return quals;
    }

    public void onTraversalDone( int[][][][] counts ) {
        out.print( header() );

        for ( int i = 0; i < 2*WIN_SIZE_SIDE; i ++ ) {
            for ( int j = 0; j < 2*WIN_SIZE_SIDE; j ++ ) {
                out.print( formatData( counts[i][j], i, j ) );
            }
        }
    }

    public Pair<Integer,Integer> countNumDeviantBases(SAMRecord read, int offset) {
        int start,end;
        if ( offset - WIN_SIZE_SIDE < 0 ) {
            start = 0;
        } else {
            start = offset - WIN_SIZE_SIDE;
        }

        if ( offset + WIN_SIZE_SIDE > read.getReadLength() ) {
            end = read.getReadLength();
        } else {
            end = offset + WIN_SIZE_SIDE;
        }

        int largeCounts = 0;
        int smallCounts = 0;
        byte[] quals = read.getBaseQualities();
        for ( int i = start; i < end; i ++ ) {
            if ( Math.abs(quals[i] - quals[offset]) >= LARGE_DEVIATION ) {
                largeCounts ++;
            } else if ( Math.abs( quals[i] - quals[offset] ) >= SMALL_DEVIATION ) {
                smallCounts ++;
            }
        }

        return new Pair<Integer,Integer>(smallCounts,largeCounts);
    }

    public boolean isMismatch( SAMRecord read, int offset, ReferenceContext ref ) {
        return ( Character.toUpperCase((char)read.getReadBases()[offset]) != ref.getBase() );
    }

    public String header() {
        return String.format(HEADER_FORMAT, "N_SMALL_DEV", "N_LARGE_DEV", "MATCH", "MISMATCH", "EXPECTED_ERR", "EMPIRICAL_ERR", "EXP_AS_Q", "EMP_AS_Q");
    }

    public String formatData( int[][] qScores, int smDev, int lgDev ) {
        int match = 0;
        int mismatch = 0;
        double expErr = 0.0;

        for ( int i = 0; i < MAX_Q_SCORE; i ++ ) {
            match += qScores[i][MATCH_OFFSET];
            mismatch += qScores[i][MM_OFFSET];
            expErr += QualityUtils.qualToProb(i)*qScores[i][MATCH_OFFSET];
            expErr += QualityUtils.qualToProb(i)*qScores[i][MM_OFFSET];
        }

        expErr = expErr/(match + mismatch);
        int expAsQ = QualityUtils.probToQual(expErr);

        double empErr = ((double)mismatch)/(match+mismatch);
        int empAsQ = QualityUtils.probToQual(empErr);

        return String.format(DATA_FORMAT, smDev, lgDev, match, mismatch, expErr, empErr, expAsQ, empAsQ);
    }
}

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
 * Time: 1:13:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class NQSClusteredZScoreWalker extends LocusWalker<LocalMapType, int[][][]> {
    static final int WIN_SIDE_SIZE = 5;
    static final int Z_SCORE_MAX = 7;
    static final int Z_SCORE_MULTIPLIER = 30; // bins are Z_SCORE * (this) rounded to the nearst int
    static final int MM_OFFSET = 1;
    static final int MATCH_OFFSET = 0;
    static final int MAX_Q_SCORE = 2 + QualityUtils.MAX_REASONABLE_Q_SCORE;

    protected int WINDOW;

    public void initialize() {
        WINDOW = 2*WIN_SIDE_SIZE+1;
    }

    public int[][][] reduceInit() {
        int[][][] q = new int[Z_SCORE_MAX*Z_SCORE_MULTIPLIER+1][MAX_Q_SCORE][2];
        for ( int i = 0; i < Z_SCORE_MAX*Z_SCORE_MULTIPLIER+1; i ++ ) {
            for ( int j = 0; j < MAX_Q_SCORE; j ++ ) {
                for ( int k = 0; k < 2; k ++ ) {
                    q[i][j][k] = 0;
                }
            }
        }

        return q;
    }


    public LocalMapType map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        return new LocalMapType(context, ref, tracker);
    }

    public int[][][] reduce(LocalMapType map, int[][][] qScores) {
        for ( int i = 0; i < map.numReads(); i ++ ) {
            SAMRecord read = map.context.getReads().get(i);
            int offset = map.context.getOffsets().get(i);
            if ( isMismatch(read, offset, map.ref) ) {
                qScores[calcZScoreBin(read,offset)][getQScoreAsInt(read,offset)][MM_OFFSET]++;
            } else {
                qScores[calcZScoreBin(read,offset)][getQScoreAsInt(read,offset)][MATCH_OFFSET]++;
            }
        }

        return qScores;
    }

    public void onTraversalDone( int[][][] zScoreBins ) {
        out.print( header() );
        for ( int i = 0; i < Z_SCORE_MAX*Z_SCORE_MULTIPLIER; i ++ ) {
            out.print( formatData(zScoreBins[i],i) );
        }
    }

    public boolean isMismatch( SAMRecord read, int offset, ReferenceContext ref ) {
        return (Character.toUpperCase((char) read.getReadBases()[offset]) != ref.getBase());
    }

    public int getQScoreAsInt( SAMRecord read, int offset ) {
        return (int) read.getBaseQualities()[offset];
    }

    public int calcZScoreBin( SAMRecord read, int offset ) {
        Pair<Double,Double> meanVar = calcWindowMeanVariance(read, offset);
        double rawZ = Math.abs((getQScoreAsInt(read, offset) - meanVar.first)/Math.sqrt(meanVar.second));
        int zBin;
        if ( rawZ < Z_SCORE_MAX ) {
            zBin = (int) Math.floor(rawZ*Z_SCORE_MULTIPLIER);
        } else {
            zBin = Z_SCORE_MULTIPLIER*Z_SCORE_MAX;
        }

        return zBin;        
    }

    public Pair<Double,Double> calcWindowMeanVariance( SAMRecord read, int offset ) {
        int start,end;
        if ( offset - WIN_SIDE_SIZE < 0 ) {
            start = 0;
        } else {
            start = offset - WIN_SIDE_SIZE;
        }

        if ( offset + WIN_SIDE_SIZE > read.getReadLength() ) {
            end = read.getReadLength();
        } else {
            end = offset + WIN_SIDE_SIZE;
        }

        double mean = 0.0;
        int n = end-start;
        byte[] quals = read.getBaseQualities();

        for ( int i = start; i < end; i ++ ) {
            if ( i != offset ) {
                mean += ((double) quals[i])/n;
            }
        }

        double var = 0.0;

        for ( int i = start; i < end; i ++ ) {
            if ( i != offset ) {
                var += Math.pow(((double) quals[i] - mean),2)/n;
            }
        }

        return new Pair<Double,Double>(mean, var);
    }

    public String header() {
        String format = "%s\t%s\t%s\t%s\t%s\t%s%n";
        return String.format(format, "ZScore", "N_obs", "Expected_Mismatch", "Empirical_Mismatch", "Expected_MM_As_Q", "Empirical_MM_As_Q");
    }

    public String formatData ( int[][] matchMismatchQ, int zScoreBin ) {
        double zScore = ( (double) zScoreBin )/Z_SCORE_MULTIPLIER;
        int match = 0;
        int mismatch = 0;
        double expMMR = 0.0;

        for ( int i = 0; i < MAX_Q_SCORE; i ++ ) {
            match += matchMismatchQ[i][MATCH_OFFSET];
            mismatch += matchMismatchQ[i][MM_OFFSET];
            expMMR += QualityUtils.qualToErrorProb((byte)i)*matchMismatchQ[i][MATCH_OFFSET];
            expMMR += QualityUtils.qualToErrorProb((byte)i)*matchMismatchQ[i][MM_OFFSET];
        }
        
        expMMR = (expMMR / ( match + mismatch ));
        double empMMR = ((double) mismatch)/(match + mismatch);
        int expMMRAsQ = QualityUtils.probToQual(1-expMMR);
        int empMMRAsQ = QualityUtils.probToQual(1-empMMR);

        String format = "%f\t%d\t%f\t%f\t%d\t%d%n";
        
        return String.format(format, zScore, match+mismatch, expMMR, empMMR, expMMRAsQ, empMMRAsQ);
    }

}

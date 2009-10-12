package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.playground.gatk.walkers.Recalibration.LocalMapType;

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
    static final int Z_SCORE_MIN = -6;
    static final int Z_SCORE_MAX = 6;
    static final int Z_SCORE_MULTIPLIER = 20; // bins are Z_SCORE * (this) rounded to the nearst int
    static final int MM_OFFSET = 1;
    static final int MATCH_OFFSET = 0;
    static final int MAX_Q_SCORE = 2 + QualityUtils.MAX_REASONABLE_Q_SCORE;

    protected int WINDOW;
    protected int Z_SCORE_RANGE;

    public void initialize() {
        WINDOW = 2*WIN_SIDE_SIZE+1;
        Z_SCORE_RANGE = Z_SCORE_MAX - Z_SCORE_MIN;
    }

    public int[][][] reduceInit() {
        int[][][] q = new int[Z_SCORE_RANGE*Z_SCORE_MULTIPLIER+1][MAX_Q_SCORE][2];
        for ( int i = 0; i < Z_SCORE_RANGE*Z_SCORE_MULTIPLIER+1; i ++ ) {
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
        for ( int i = 0; i < Z_SCORE_RANGE*Z_SCORE_MULTIPLIER; i ++ ) {
            for ( int j = 0; j < MAX_Q_SCORE; j ++ ) {
                out.print( formatData(zScoreBins[i][j], i, j) );
            }
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
        double rawZ = (getQScoreAsInt(read, offset) - meanVar.first)/Math.sqrt(meanVar.second);
        int zBin;
        if ( rawZ >= Z_SCORE_MAX ) {
            zBin = (int) Math.floor(Z_SCORE_RANGE*Z_SCORE_MULTIPLIER);
        } else if ( rawZ <= Z_SCORE_MIN ) {
            zBin = 0;
        } else if ( rawZ > 0 ) {
            zBin = (int) Math.floor(rawZ*Z_SCORE_MULTIPLIER) - (int) Math.floor(Z_SCORE_MIN*Z_SCORE_MULTIPLIER);
        } else {
            zBin = (int) Math.floor(-Z_SCORE_MIN*Z_SCORE_MULTIPLIER) + (int) Math.floor(rawZ*Z_SCORE_MULTIPLIER);
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

    public String formatData ( int[] matchMismatch, int zScoreBin, int q ) {
        String format = "%f\t%d\t%f\t%f\t%d\t%d%n";

        double zScore = ( (double) zScoreBin )/Z_SCORE_MULTIPLIER;
        int counts = matchMismatch[MATCH_OFFSET] + matchMismatch[MM_OFFSET];
        double empMMR = (((double)matchMismatch[MM_OFFSET])/counts);
        double expMMR = QualityUtils.qualToErrorProb((byte) q);
        int empMMRAsQ = QualityUtils.probToQual(1-empMMR);
        
        return String.format(format, zScore, counts, expMMR, empMMR, q, empMMRAsQ);
    }

}

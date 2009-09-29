package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Pair;

import java.util.List;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 28, 2009
 * Time: 12:02:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class QualityScoreDistributionWalker extends LocusWalker<List<Pair<Integer,Integer>>, int[][]> {
    @Argument(fullName="maxReadLength", shortName="rl", doc="Maximum read length in the bam file", required=true)
    protected int MAX_READ_LENGTH = 125;


    protected final int MAX_REASONABLE_Q_SCORE = 3 + QualityUtils.MAX_REASONABLE_Q_SCORE;
    public void initialize() {

    }

    public int[][] reduceInit() {
        int [][] qualityCounts = new int[MAX_READ_LENGTH][MAX_REASONABLE_Q_SCORE];
        for ( int i = 0; i < MAX_READ_LENGTH; i ++ ) {
            for ( int j = 0; j < MAX_REASONABLE_Q_SCORE; j ++ ) {
                qualityCounts[i][j] = 0;
            }
        }

        return qualityCounts;
    }

    public List<Pair<Integer,Integer>> map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        LinkedList<Pair<Integer,Integer>> qualByPos = new LinkedList<Pair<Integer,Integer>>();
        for ( int r = 0; r < context.numReads(); r ++ ) {
            if( context.getReads().get(r).getReadNegativeStrandFlag() ) {
                int of = context.getReads().get(r).getReadLength() - context.getOffsets().get(r);
                int qu = context.getReads().get(r).getBaseQualities()[context.getOffsets().get(r)];
                qualByPos.add(new Pair<Integer,Integer>(of, qu));
            }
            else {
                int of = context.getOffsets().get(r);
                int qu = context.getReads().get(r).getBaseQualities()[of];
                qualByPos.add(new Pair<Integer,Integer>(of, qu));
            }
        }

        return qualByPos;
    }

    public int[][] reduce( List<Pair<Integer,Integer>> qualByPos, int[][] counts ) {
        for ( Pair<Integer,Integer> qPosPair : qualByPos ) {
            counts[qPosPair.getFirst()][qPosPair.getSecond()] ++;
        }

        return counts;
    }

    public void onTraversalDone( int[][] counts ) {
        for ( int pos = 0; pos < MAX_READ_LENGTH; pos ++ ) {
            String outStr = Integer.toString(pos);
            for ( int qual = 0; qual < MAX_REASONABLE_Q_SCORE; qual ++) {
                outStr = outStr + "\t" + Integer.toString(counts[pos][qual]);
            }
            out.printf("%s%n", outStr);
        }
    }
}

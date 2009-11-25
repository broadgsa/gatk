package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 20, 2009
 * Time: 9:41:57 AM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value={DataSource.REFERENCE}, referenceMetaData = {@RMD(name="dbsnp",type= rodDbSNP.class)})
public class MinimumNQSWalker extends LocusWalker<Pair<List<Pair<Integer,Integer>>,List<Pair<Integer,Integer>>>, int[][][]> {
        private static int MM_OFFSET = 1;
    private static int MATCH_OFFSET = 0;
    private static int QSCORE_MAX = 1 + QualityUtils.MAX_REASONABLE_Q_SCORE;
    @Argument(fullName="windowSize", shortName="ws", doc="Size of the window (in one direction)", required=true)
    private int winSide = 4;

    public void initialize() {
        out.printf("%s%n", makeHeader());
    }

    public int[][][] reduceInit() {
        int[][][] counts = new int[QSCORE_MAX][QSCORE_MAX][2];
        for ( int i = 0; i < QSCORE_MAX; i ++ ) {
            for ( int j = 0; j < QSCORE_MAX; j ++ ) {
                counts[i][j][1]=0;
                counts[i][j][0]=0;
            }
        }
        return counts;
    }

    public int[][][] reduce(Pair<List<Pair<Integer,Integer>>,List<Pair<Integer,Integer>>> map, int[][][] prevReduce) {
        if ( map != null ) {
            List<Pair<Integer,Integer>> matchingQualityNQSPairs = map.getFirst();
            List<Pair<Integer,Integer>> mismatchingQualityNQSPairs = map.getSecond();
            for ( Pair<Integer,Integer> p : matchingQualityNQSPairs ) {
                prevReduce[p.getFirst()][p.getSecond()][MATCH_OFFSET] ++;
            }

            for ( Pair<Integer,Integer> p : mismatchingQualityNQSPairs ) {
                prevReduce[p.getFirst()][p.getSecond()][MM_OFFSET] ++;
            }
        }

        return prevReduce;
    }

    public Pair<List<Pair<Integer,Integer>>,List<Pair<Integer,Integer>>> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ArrayList<Pair<Integer,Integer>> matchingQualityNQSPairs = new ArrayList<Pair<Integer,Integer>>();
        ArrayList<Pair<Integer,Integer>> mismatchingQualityNQSPairs = new ArrayList<Pair<Integer,Integer>>();
        if ( (Variation) tracker.lookup("dbsnp",null) == null ) {
            for ( int r = 0; r < context.size(); r ++ ) {
                SAMRecord read = context.getReads().get(r);
                int offset = context.getOffsets().get(r);
                int quality = read.getBaseQualities()[offset];
                int NQS = getNghdQ(offset, read);
                Pair<Integer,Integer> qualityNQSPair = new Pair<Integer,Integer> (quality,NQS);
                if ( BaseUtils.basesAreEqual(read.getReadBases()[offset], (byte) ref.getBase()) ) {
                    matchingQualityNQSPairs.add(qualityNQSPair);
                } else {
                    mismatchingQualityNQSPairs.add(qualityNQSPair);
                }
            }

            return new Pair<List<Pair<Integer,Integer>>,List<Pair<Integer,Integer>>>(matchingQualityNQSPairs,mismatchingQualityNQSPairs);
        } else {
            return null;
        }
    }


    public void onTraversalDone( int[][][] reduce ) {
        for ( int qc = 0; qc < QSCORE_MAX; qc ++ ) {
            for ( int qn = 0; qn < QSCORE_MAX; qn ++ ) {
                out.printf("%s%n", formatData(reduce[qc][qn],qc,qn));
            }
        }
    }

    public int getNghdQ(int off, SAMRecord read) {
        // System.out.println("getNghdQ");
        byte minQ = Byte.MAX_VALUE;
        byte[] quals = read.getBaseQualities();
        int rdlnth = read.getReadLength();
        int start;
        int end;
        if ( off - winSide < 0 ) {
            start = 0;
        } else {
            start = off - winSide;
        }

        if ( off + winSide > rdlnth ) {
            end = rdlnth;
        } else {
            end = off + winSide;
        }

        for ( int i = start; i < end; i ++ ) {
            if ( i != off ) {
                byte q = quals[i];
                if ( q < minQ ) {
                    minQ = q;
                }
            }
        }

        return minQ;
    }

    private String makeHeader() {
        return String.format("%s\t%s\t%s\t%s\t%s\t%s","Reported_Q","Min_Nghd_Q","N_observations","Mm_rate","Empirical_Q","Q_diff");
    }

    private String formatData( int[] mmArray, int qCenter, int qNghd ) {
        int counts = mmArray[MM_OFFSET]+mmArray[MATCH_OFFSET];
        double mismatch = ((double)mmArray[MM_OFFSET]/counts);
        byte qEmp = QualityUtils.probToQual(1-mismatch);
        return String.format("%d\t%d\t%d\t%f\t%d\t%d", qCenter, qNghd, counts, mismatch, qEmp,qEmp-qCenter);
    }
}


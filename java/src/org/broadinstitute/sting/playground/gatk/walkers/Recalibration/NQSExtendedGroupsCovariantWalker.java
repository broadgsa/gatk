package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.List;

import net.sf.samtools.SAMRecord;


/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 23, 2009
 * Time: 11:26:52 AM
 * To change this template use File | Settings | File Templates.
 */
public class NQSExtendedGroupsCovariantWalker extends LocusWalker<LocalMapType, long[][][]> {

    final int NEIGHBORHOOD_SIZE = 5;
    final int MM_OFFSET = 1;
    final int COUNT_OFFSET = 0;
    final int MAX_QSCORE = QualityUtils.MAX_REASONABLE_Q_SCORE+1;
    int NQS_GROUPS;

    public static final String DATA_FORMAT = "%d\t%d\t%d\t%d\t%d%n";
    public static final String TEXT_FORMAT = "%s\t%s\t%s\t%s\t%s%n";


    // Walker methods

    public void initialize() {
        NQS_GROUPS = (MAX_QSCORE)*(MAX_QSCORE+1)/2 + MAX_QSCORE;
        // this is from solving the recurrent equation
        // a_{n+1} = a_n + n + 1
    }

    public long[][][] reduceInit() {
        long[][][] covariantCounts = new long[NQS_GROUPS][2][MAX_QSCORE];

        for ( int i = 0; i < NQS_GROUPS; i ++ ) {
            for ( int j = 0; j < 2; j++ ) {
                for ( int k = 0; k < MAX_QSCORE; k ++) {
                    covariantCounts[i][j][k] = 0;
                }
            }
        }

        return covariantCounts;
    }

    public LocalMapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return new LocalMapType(context, ref, tracker);
    }

    public long[][][] reduce( LocalMapType map, long[][][] cumulativeBins ) {

        for ( int i = 0; i < map.numReads(); i ++ ) {
            int groupNumber = calcGroupNumber(i, map.context);
            int qualityAsInt = getBaseQualityAsInt(map.context,i);
            cumulativeBins[groupNumber][COUNT_OFFSET][qualityAsInt] ++;
            if ( ! isRef(i, map.context, map.ref ) ) {
                cumulativeBins[groupNumber][MM_OFFSET][qualityAsInt] ++;
            }
        }

        return cumulativeBins;
    }

    public void onTraversalDone(long[][][] cumulativeBins) {
        printNQSMismatchTable(cumulativeBins);
    }

    // Map-Reduce helper methods

    public int calcGroupNumber( int read, AlignmentContext context ) {
        return calcGroupNumberFromQualities(getBaseQualityAsInt(context, read), getMinNeighborhoodBaseQualityAsInt(context, read));
    }

    /*
     * calculates group number as
     * #:   nghd left       base         nghd right
     * 0:     000             0             000
     * 1:     000             1             000
     * 2:     111             1             111
     * 3:     000             2             000
     * 4:     111             2             111
     * 5:     222             2             222
     * 6:     000             3             000
     * 7:     111             3             111
     * etc
     */
    public int calcGroupNumberFromQualities( int baseQuality, int minNeighborhoodQuality ) {
        int groupNumber = (baseQuality)*(baseQuality+1)/2;
        if ( minNeighborhoodQuality <= baseQuality ) {
            groupNumber += minNeighborhoodQuality;
        } else {
            groupNumber += baseQuality; // yes we want to do this
        }

        return groupNumber;
    }

    public int getBaseQualityAsInt( AlignmentContext context, int read ) {
        return (int) context.getReads().get(read).getBaseQualities()[context.getOffsets().get(read)];
    }

    public int getMinNeighborhoodBaseQualityAsInt( AlignmentContext context, int read ) {
        byte[] qualities = context.getReads().get(read).getBaseQualities();
        int offset = context.getOffsets().get(read);
        int readLength = qualities.length;
        int start = (offset - NEIGHBORHOOD_SIZE < readLength) ? 0 : offset - NEIGHBORHOOD_SIZE;
        int end = (offset + NEIGHBORHOOD_SIZE > readLength) ? readLength : offset+NEIGHBORHOOD_SIZE;
        byte minQuality = Byte.MAX_VALUE;
        
        for ( int i = start; i < end; i ++ ) {
            if ( i != offset ) {
                if ( qualities[i] < minQuality ) {
                    minQuality = qualities[i];
                }
            }
        }

        return (int) minQuality;
    }

    public boolean isRef(int read, AlignmentContext context, ReferenceContext ref ) {
        return (context.getReads().get(read).getReadBases()[context.getOffsets().get(read)] == ref.getBase());
    }

    // printing helper methods

    public void printNQSMismatchTable(long[][][] cumulativeBins) {
        out.print(createHeader());
        for ( int qscore = 0; qscore < MAX_QSCORE; qscore ++ ) {
            for ( int group = 0; group < NQS_GROUPS; group ++ ) {
                out.print(formatNQSMismatchCountString(qscore,group,cumulativeBins));
            }
        }
    }

    public String createHeader() {
        return String.format(TEXT_FORMAT, "Qscore_Reported", "NQS_Group", "Mismatches", "Total", "Empirical_Qscore");
    }

    public String formatNQSMismatchCountString(int qscore, int group, long[][][] cumulativeBins) {
        long mm = cumulativeBins[group][MM_OFFSET][qscore];
        long ct = cumulativeBins[group][COUNT_OFFSET][qscore];
        double mmr = (ct > 0) ? mm/ct : -1;

        return String.format(DATA_FORMAT, qscore, group, mm, ct, QualityUtils.probToQual(1-mmr));
    }
}

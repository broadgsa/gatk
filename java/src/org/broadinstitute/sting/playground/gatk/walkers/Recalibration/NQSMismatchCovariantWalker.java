package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.playground.gatk.walkers.Recalibration.LocalMapType;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.QualityUtils;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 14, 2009
 * Time: 1:37:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class NQSMismatchCovariantWalker extends LocusWalker<LocalMapType, long[][][]> {
    public static final int NQS_GROUPS = 35;
    public static final int NQS_RESOLUTION = 1;
    public static final int NQS_DIFFERENCE = 5;
    public static final int NEIGHBORHOOD_SIZE = 5;
    public static final int NQS_QSCORE_START = 5;
    public static final int MM_OFFSET = 0;
    public static final int COUNT_OFFSET = 1;
    public static final int FAIL_OFFSET = 2;
    public static final String DATA_FORMAT = "%d                       %d                             %d                               %d                                 %d                            %d%n";
    public static final String HEADER_FORMAT = "%s    %s     %s     %s     %s     %s%n";

    public long[][][] reduceInit() {
        long[][][] mismatchesByNQS = new long[(int) QualityUtils.MAX_QUAL_SCORE+1][this.numNQSGroups()+1][2];
        for ( int qualityScore = 0; qualityScore <= (int) QualityUtils.MAX_QUAL_SCORE; qualityScore ++) {
            for ( int nqsGroup = 0; nqsGroup <= numNQSGroups(); nqsGroup++ ) {
                mismatchesByNQS[qualityScore][nqsGroup][MM_OFFSET] = 0;
                mismatchesByNQS[qualityScore][nqsGroup][COUNT_OFFSET] = 0;
            }
        }
        return mismatchesByNQS;
    }

    public LocalMapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        return new LocalMapType(context,ref,tracker);
    }

    public long[][][] reduce(LocalMapType rawMapData, long[][][] cumMismatchArray) {
        if ( ! isDBSNP(rawMapData.tracker) ) {
            for( int read = 0; read < rawMapData.numReads(); read ++ ) {
                int bestNQSRelativeToStart = getBestNQSRelativeToStart(rawMapData.context, read, nqsDifference(),
                                               nqsResolution(), centerBaseQualityStart(), windowSize());
                if( isValidRelativeNQS(bestNQSRelativeToStart) ) {
                    cumMismatchArray[rawMapData.qScore(read)][bestNQSRelativeToStart][COUNT_OFFSET]++;
                    if( isMismatch(read, rawMapData.context, rawMapData.ref) ) {
                        cumMismatchArray[rawMapData.qScore(read)][bestNQSRelativeToStart][MM_OFFSET]++;
                    }
                }
            }
        }

        return cumMismatchArray;
    }

    public void onTraversalDone(long[][][] cumulativeCounts) {
        printNQSMismatchTable(cumulativeCounts);
    }

    public int getBestNQSRelativeToStart(AlignmentContext context, int readNum, int ctrNghdNQSDif, int nqsStepSize,
                                         int ctrBaseQStart, int windSize) {
        int minNeighborhoodQualityScore = findMinQScoreInWindow(context, readNum, windSize);
        int centerQScore = getQualityScore(context,readNum);
        int stepsToCutoffNeighborhood = (int) Math.floor((minNeighborhoodQualityScore - (ctrBaseQStart - ctrNghdNQSDif))/nqsStepSize);
        int stepsToCutoffCenter = (int) Math.floor((centerQScore-ctrBaseQStart )/nqsStepSize);

        if(stepsToCutoffNeighborhood < stepsToCutoffCenter) {
            return stepsToCutoffNeighborhood;
        } else {
            return stepsToCutoffCenter;
        }
    }

    protected void printNQSMismatchTable(long[][][] cumulativeCounts) {
        out.printf("%s", createHeader() );
        for ( int qscore = 0; qscore <= QualityUtils.MAX_QUAL_SCORE; qscore ++ ) {
        for( int nqsGroup = 0; nqsGroup <= numNQSGroups(); nqsGroup ++ ) {
            out.printf("%s", formatNQSMismatchCountString(cumulativeCounts, nqsGroup, qscore));
        }
        }
    }

    protected String createHeader() {
        return String.format(HEADER_FORMAT,"Qscore Reported","Qscore Threshold at Locus", "Minimum Neighborhood Quality", "NQS Group",
                             "# Non-dbSNP Mismatches", "Total Non-dbSNP Counts");
    }

    protected String formatNQSMismatchCountString(long[][][] counts, int nqsGroup, int qscore) {
        return String.format(DATA_FORMAT,qscore,getCenterThreshold(nqsGroup),getNeighborhoodThreshold(nqsGroup),nqsGroup,
                             counts[qscore][nqsGroup][MM_OFFSET], counts[qscore][nqsGroup][COUNT_OFFSET]);
    }

    protected int getCenterThreshold(int nqsGroup) {
        return centerBaseQualityStart() + nqsGroup*nqsResolution();
    }

    protected int getNeighborhoodThreshold(int nqsGroup) {
        return getCenterThreshold(nqsGroup)-nqsDifference();
    }

    protected int numNQSGroups() {
        return NQS_GROUPS;
    }

    protected int nqsResolution() {
        return NQS_RESOLUTION;
    }

    protected int nqsDifference() {
        return NQS_DIFFERENCE;
    }

    protected boolean isValidRelativeNQS(int relNQS) {
        return relNQS >= 0;
    }

    protected boolean isDBSNP(RefMetaDataTracker tracker) {

        return false;
        // return ( tracker.lookup("dbSNP",null) != null );
    }

    protected boolean isMismatch(int readNum, AlignmentContext context, ReferenceContext ref) {
        return ( (char) getBases(context,readNum)[ getOffset(context,readNum) ] ) != ref.getBase();
    }

    protected byte[] getBases(AlignmentContext context, int offset) {
        return context.getReads().get(offset).getReadBases();
    }

    protected int getOffset(AlignmentContext context, int read) {
        return context.getOffsets().get(read);
    }

    protected int windowSize() {
        return NEIGHBORHOOD_SIZE;
    }

    protected int centerBaseQualityStart() {
        return NQS_QSCORE_START;
    }

    protected int getQualityScore(AlignmentContext context, int readNo) {
        return getQualityScore(context.getReads().get(readNo), context.getOffsets().get(readNo));
    }

    protected int getQualityScore(SAMRecord read, int offset) {
        return (int) read.getBaseQualities()[offset];
    }

    protected int findMinQScoreInWindow(AlignmentContext context, int readNo, int windowSize) {
        SAMRecord read = context.getReads().get(readNo);
        int offset = context.getOffsets().get(readNo);
        int start;
        int end;

        if ( offset - windowSize < 0) {
            start = 0;
        } else {
            start = offset - windowSize;
        }

        if ( offset + windowSize > read.getReadLength() ) {
            end = read.getReadLength();
        } else {
            end = offset + windowSize;
        }

        int minQualityScore = Integer.MAX_VALUE;
        for (int positionInRead = start; positionInRead < end; positionInRead++ ) {
            int qScoreAtPosition = getQualityScore(read,positionInRead);
            if ( qScoreAtPosition < minQualityScore ) {
                minQualityScore = qScoreAtPosition;
            }
        }

        return minQualityScore;
    }

}


package edu.mit.broad.picard.util;

import edu.mit.broad.sam.SAMRecord;

/**
 * Utility mthods for pairs of SAMRecords
 */
public class SamPairUtil {

    // TODO: KT and TF say this is more complicated than what I have here
    public static boolean isProperPair(final SAMRecord firstEnd, final SAMRecord secondEnd, boolean jumpingLibrary) {
        if (firstEnd.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
            return false;
        }
        if (!firstEnd.getReferenceName().equals(secondEnd.getReferenceName())) {
            return false;
        }
        if (firstEnd.getReadNegativeStrandFlag() == secondEnd.getReadNegativeStrandFlag()) {
            return false;
        }
        final SAMRecord positiveEnd;
        final SAMRecord negativeEnd;
        if (firstEnd.getReadNegativeStrandFlag()) {
            positiveEnd = secondEnd;
            negativeEnd = firstEnd;
        } else {
            positiveEnd = firstEnd;
            negativeEnd = secondEnd;
        }
        if (!jumpingLibrary) {
            return positiveEnd.getAlignmentStart() < negativeEnd.getAlignmentStart() + negativeEnd.getReadBases().length;
        } else {
            return negativeEnd.getAlignmentStart() < positiveEnd.getAlignmentStart() + positiveEnd.getReadBases().length;
        }
    }

    public static int computeInsertSize(final SAMRecord firstEnd, final SAMRecord secondEnd) {
        if (firstEnd.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
            return 0;
        }
        if (!firstEnd.getReferenceName().equals(secondEnd.getReferenceName())) {
            return 0;
        }
        int firstEnd5PrimePosition = firstEnd.getReadNegativeStrandFlag()? firstEnd.getAlignmentEnd(): firstEnd.getAlignmentStart();
        int secondEnd5PrimePosition = secondEnd.getReadNegativeStrandFlag()? secondEnd.getAlignmentEnd(): secondEnd.getAlignmentStart();
        return secondEnd5PrimePosition - firstEnd5PrimePosition;
    }

    /**
     * Write the mate info for two SAMRecords
     */
    public static void setMateInfo(final SAMRecord samRecord, final SAMRecord mate) {
        if (!samRecord.getMateUnmappedFlag()) {
            samRecord.setMateReferenceName(mate.getReferenceName());
            samRecord.setMateAlignmentStart(mate.getAlignmentStart());
            samRecord.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());
        } else {
            samRecord.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            samRecord.setMateUnmappedFlag(true);
        }
        if (!mate.getMateUnmappedFlag()) {
            mate.setMateReferenceName(samRecord.getReferenceName());
            mate.setMateAlignmentStart(samRecord.getAlignmentStart());
            mate.setMateNegativeStrandFlag(samRecord.getReadNegativeStrandFlag());
        } else {
            mate.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            mate.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            mate.setMateUnmappedFlag(true);
        }
    }


}

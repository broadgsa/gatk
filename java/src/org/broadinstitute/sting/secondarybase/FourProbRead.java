package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.ArrayList;

/**
 * FourProbRead contains the four-prob information for each base in a read.
 */
public class FourProbRead extends ArrayList<FourProb> {
    /**
     * Initialize the container with the specified capacity.
     *
     * @param initialCapacity  the number of bases in the read
     */
    public FourProbRead(int initialCapacity) {
        super(initialCapacity);
    }

    /**
     * Returns a subset of the FourProbRead.
     *
     * @param cycleStart  the starting cycle (0-based, inclusive)
     * @param cycleStop   the ending cycle (0-based, inclusive)
     * @return a FourProbRead that is a subset of this FourProbRead
     */
    public FourProbRead getSubset(int cycleStart, int cycleStop) {
        FourProbRead subset = new FourProbRead(cycleStop - cycleStart + 1);

        for (int cycle = cycleStart, offset = 0; cycle <= cycleStop; cycle++, offset++) {
            subset.add(offset, this.get(cycle));
        }

        return subset;
    }

    /**
     * Get the read sequence at a specified rank.
     *
     * @param rank  the rank of the sequence to return (0=best, 1=second-best, 2=third-best, 3=fourth-best)
     * @return the read sequence at the specified rank
     */
    public String getBaseSequenceAtGivenRank(int rank) {
        String pseq = "";

        for ( FourProb fp : this ) {
            pseq += fp.baseAtRank(rank);
        }

        return pseq;
    }

    /**
     * Get the primary read sequence.
     *
     * @return  the primary read sequence
     */
    public String getPrimaryBaseSequence() {
        return getBaseSequenceAtGivenRank(0);
    }

    /**
     * Get the secondary read sequence.
     *
     * @return  the secondary read sequence
     */
    public String getSecondaryBaseSequence() {
        return getBaseSequenceAtGivenRank(1);
    }

    /**
     * Get the SAM spec-conformant SQ tag that will be written to the SAM/BAM file.
     *
     * @param rr  the raw read
     * @return the byte array for the SQ tag (first two bits: the base identity, the last six bits: -10*log10(p3/p2)
     */
    public byte[] getSQTag(RawRead rr) {
        byte[] sqtag = new byte[this.size()];

        for (int cycle = 0; cycle < this.size(); cycle++) {
            FourProb fp = this.get(cycle);

            int fpPrimaryBaseIndex = fp.indexAtRank(0);
            int rawPrimaryBaseIndex = BaseUtils.simpleBaseToBaseIndex(rr.getSequenceAsString().charAt(cycle));

            int fpSecondaryBaseIndex = (fpPrimaryBaseIndex == rawPrimaryBaseIndex) ? fp.indexAtRank(1) : fpPrimaryBaseIndex;

            double qualdiff = -10.0*Math.log10(fp.probAtRank(2)/fp.probAtRank(1));

            sqtag[cycle] = QualityUtils.baseAndProbDiffToCompressedQuality(fpSecondaryBaseIndex, qualdiff);
        }

        return sqtag;
    }
}

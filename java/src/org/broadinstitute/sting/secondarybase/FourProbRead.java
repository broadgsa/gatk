package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.ArrayList;

/**
 * FourProbRead contains the four-prob information for each base in a read.
 */
public class FourProbRead extends ArrayList<FourProb> {
    /**
     * Initialize the container with the specified capacity
     *
     * @param initialCapacity  the number of bases in the read
     */
    public FourProbRead(int initialCapacity) {
        super(initialCapacity);
    }

    /**
     * Get the read sequence at a specified rank
     *
     * @param rank  the rank of the sequence to return (0=best, 1=second-best, 2=third-best, 3=fourth-best)
     * @return the read sequence at the specified rank
     */
    public String getBaseSequenceAtGivenRank(int rank) {
        String pseq = "";

        for (int cycle = 0; cycle < this.size(); cycle++) {
            FourProb fp = this.get(cycle);

            pseq += fp.baseAtRank(rank);
        }

        return pseq;
    }

    /**
     * Get the primary read sequence
     * @return  the primary read sequence
     */
    public String getPrimaryBaseSequence() {
        return getBaseSequenceAtGivenRank(0);
    }

    /**
     * Get the secondary read sequence
     * @return  the secondary read sequence
     */
    public String getSecondaryBaseSequence() {
        return getBaseSequenceAtGivenRank(1);
    }

    /**
     * Get the SAM spec-conformant SQ tag that will be written to the SAM/BAM file.
     *
     * @param rr  the raw read
     * @return the byte array for the SQ tag (first two bits: the base identity, the last six bits: -10*log10(p3/p2).
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

            //System.out.println("SQTAG: " + qualdiff + " " + sqtag[cycle] + " " + QualityUtils.compressedQualityToProbDiff(sqtag[cycle]));
        }

        return sqtag;
    }
}

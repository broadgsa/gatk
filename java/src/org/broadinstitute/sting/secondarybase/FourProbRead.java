package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.ArrayList;

public class FourProbRead extends ArrayList<FourProb> {
    public FourProbRead(int initialCapacity) {
        super(initialCapacity);
    }

    public String getBaseSequenceAtGivenRank(int rank) {
        String pseq = "";

        for (int cycle = 0; cycle < this.size(); cycle++) {
            FourProb fp = this.get(cycle);

            pseq += fp.baseAtRank(rank);
        }

        return pseq;
    }

    public String getPrimaryBaseSequence() {
        return getBaseSequenceAtGivenRank(0);
    }

    public String getSecondaryBaseSequence() {
        return getBaseSequenceAtGivenRank(1);
    }

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

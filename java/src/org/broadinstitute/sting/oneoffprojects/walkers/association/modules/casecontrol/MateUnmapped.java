package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/2/11
 * Time: 2:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateUnmapped extends ProportionTest {
    private final int MAPQ_THRESHOLD = 5;

    public Pair<Number,Number> map(ReadBackedPileup pileup) {
        int numMatedReads = 0;
        int numPairUnmapped = 0;
        for (PileupElement e : pileup ) {
            if (e.getRead().getReadPairedFlag() && e.getMappingQual() >= MAPQ_THRESHOLD ) {
                ++numMatedReads;
                if ( e.getRead().getMateUnmappedFlag() ) {
                    ++numPairUnmapped;
                }
            }
        }

        return new Pair<Number,Number>(numPairUnmapped,numMatedReads);
    }

    public boolean usePreviouslySeenReads() { return false; }
}

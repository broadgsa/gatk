package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProperPairs extends ZStatistic {

    public Pair<Number,Number> map(ReadBackedPileup rbp) {
        int numReads = 0;
        int numPropPair = 0;
        for (PileupElement e : rbp ) {
            ++numReads;
            if ( e.getRead().getReadPairedFlag() && e.getRead().getProperPairFlag() ) {
                ++numPropPair;
            }
        }

        return new Pair<Number,Number>(numPropPair,numReads);
    }

    public boolean usePreviouslySeenReads() { return false; }
}

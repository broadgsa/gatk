package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.ZStatistic;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateSameStrand extends ZStatistic {

    public Pair<Number,Number> map(ReadBackedPileup rbp) {
        int numPairs = 0;
        int mateSameStrand = 0;
        for (PileupElement e : rbp ) {
            if ( e.getRead().getReadPairedFlag() ) {
                ++numPairs;
                if ( e.getRead().getMateNegativeStrandFlag() == e.getRead().getReadNegativeStrandFlag() ) {
                    ++mateSameStrand;
                }
            }
        }

        return new Pair<Number,Number>(mateSameStrand,numPairs);
    }

    public int getWindowSize() { return 100; }
    public int slideByValue() { return 10; }
    public boolean usePreviouslySeenReads() { return false; }
}

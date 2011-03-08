package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.UStatistic;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateMappingQuality extends UStatistic {

    public Collection<Number> map(ReadBackedPileup rbp) {
        ArrayList<Integer> mateMapQ = new ArrayList<Integer>(rbp.size());
        for ( PileupElement e : rbp ) {
            if ( e.getRead().getReadPairedFlag() && e.getRead().getAttribute("MQ") != null) {
                mateMapQ.add((Integer) e.getRead().getAttribute("MQ"));
            }
        }

        return (Collection) mateMapQ;
    }

    public int getWindowSize() { return 100; }
    public int slideByValue() { return 10; }
    public boolean usePreviouslySeenReads() { return false; }
}

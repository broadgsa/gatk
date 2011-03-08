package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.TStatistic;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.UStatistic;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class BaseQualityScore extends TStatistic {

    public Collection<Number> map(ReadBackedPileup rbp) {
        ArrayList<Integer> baseQuals = new ArrayList<Integer>(rbp.size());
        for (PileupElement e : rbp ) {
            baseQuals.add((int)e.getQual());
        }

        return (Collection) baseQuals;
    }

    public int getWindowSize() { return 100; }
    public int slideByValue() { return 10; }
    public boolean usePreviouslySeenReads() { return true; }
}

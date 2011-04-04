package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class MappingQuality0 extends ZStatistic {

    public Pair<Number,Number> map(ReadBackedPileup rbp) {
        int total = 0;
        int mq0 = 0;
        for ( PileupElement e : rbp ) {
            ++total;
            if (e.getMappingQual() == 0)
                ++mq0;
        }

        return new Pair<Number,Number>(mq0,total);
    }

    public boolean usePreviouslySeenReads() { return false; }

}

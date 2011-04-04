package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/1/11
 * Time: 3:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateOtherContig extends ZStatistic {

    public boolean usePreviouslySeenReads() { return false; }

    public Pair<Number,Number> map(ReadBackedPileup pileup) {
        int tot = 0;
        int otherCon = 0;
        for ( PileupElement e : pileup ) {
            if ( e.getRead().getReadPairedFlag() ) {
                ++tot;
                if ( ! e.getRead().getMateReferenceIndex().equals(e.getRead().getReferenceIndex()) ) {
                    ++otherCon;
                }
            }
        }

        return new Pair<Number,Number>(otherCon,tot);
    }

}

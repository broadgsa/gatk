package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import ch.qos.logback.classic.filter.ThresholdFilter;
import org.broadinstitute.sting.oneoffprojects.walkers.association.RegionalAssociationWalker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: 4/9/11
 * Time: 3:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadsAberrantInsertSize extends ProportionTest {

    private int HI_THRESHOLD;
    private int LO_THRESHOLD;

    @Override
    public Pair<Number,Number> map(ReadBackedPileup rbp) {
        int total = 0;
        int wonky = 0;
        for (PileupElement e : rbp ) {
            if ( e.getRead().getReadPairedFlag() ) {
                ++total;
                if ( e.getRead().getInferredInsertSize() > HI_THRESHOLD || e.getRead().getInferredInsertSize() < LO_THRESHOLD ) {
                    ++wonky;
                }
            }
        }

        return new Pair<Number,Number>(wonky,total);
    }

    public void init(RegionalAssociationWalker parent) {
        super.init(parent);
        HI_THRESHOLD = 200;
	LO_THRESHOLD = 50;
    }

    public boolean usePreviouslySeenReads() { return false; }
}

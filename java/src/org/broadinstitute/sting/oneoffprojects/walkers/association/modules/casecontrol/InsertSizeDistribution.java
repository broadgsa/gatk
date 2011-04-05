package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author chartl
 */
public class InsertSizeDistribution extends ValueTest {

    public boolean usePreviouslySeenReads() { return false; }

    public Collection<Number> map(ReadBackedPileup pileup) {
        List<Integer> insertSizes = new ArrayList<Integer>(pileup.size());
        for ( PileupElement e : pileup ) {
            insertSizes.add(e.getRead().getInferredInsertSize());
        }

        return (Collection) insertSizes;
    }

}

package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadMappingQuality extends TStatistic {

    public Collection<Number> map(ReadBackedPileup rbp) {
        ArrayList<Integer> mapQuals = new ArrayList<Integer>(rbp.size());
        for ( PileupElement e : rbp ) {
            mapQuals.add(e.getMappingQual());
        }

        return (Collection) mapQuals;
    }

    public boolean usePreviouslySeenReads() { return false; }
}

package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.MapExtender;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/8/11
 * Time: 1:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceMismatches extends ZStatistic {

    final static int[] BASE_INDEX = {0,1,2,3};

    int currentRefBase = 0;

    @Override
    public Map<Sample,Object> mapLocus(MapExtender extender) {
        currentRefBase = BaseUtils.simpleBaseToBaseIndex(extender.getReferenceContext().getBase());
        return super.mapLocus(extender);
    }

    public Pair<Number,Number> map(ReadBackedPileup rbp) {
        int[] counts = rbp.getBaseCounts();
        int total = 0;
        int nonref = 0;
        for ( int base : BASE_INDEX ) {
            total += counts[base];
            if ( base != currentRefBase ) {
                nonref += counts[base];
            }
        }

        return new Pair<Number,Number>(nonref,total);
    }

    public boolean usePreviouslySeenReads() { return true; }
}

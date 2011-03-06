package org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 12:58:16 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TStatistic extends CaseControl<Collection<Number>> {

    public abstract Collection<Number> map(ReadBackedPileup rbp );

    public Collection<Number> add(Collection<Number> left, Collection<Number> right) {
        if ( left instanceof ArrayList) {
            ((ArrayList) left).addAll(right);
            return left;
        } else if ( left instanceof Set ) {
            ((Set) left).addAll(right);
            return left;
        } else {
            List<Number> newList = new ArrayList<Number>(left.size()+right.size());
            newList.addAll(left);
            newList.addAll(right);
            return newList;
        }
    }

}
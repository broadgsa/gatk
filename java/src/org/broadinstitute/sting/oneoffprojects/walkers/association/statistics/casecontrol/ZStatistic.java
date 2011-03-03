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
 * Time: 12:58:38 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ZStatistic extends CaseControl<Pair<Number,Number>> {
    public abstract Pair<Number,Number> map(ReadBackedPileup rbp );

    public Pair<Number,Number> add(Pair<Number,Number> left, Pair<Number,Number> right) {
        return new Pair<Number,Number>(left.first.doubleValue()+right.first.doubleValue(),
                                       left.second.doubleValue()+right.second.doubleValue());
    }

}
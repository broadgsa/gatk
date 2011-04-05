package org.broadinstitute.sting.oneoffprojects.walkers.association.statistics;

import org.broadinstitute.sting.utils.collections.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: 4/5/11
 * Time: 12:26 PM
 * To change this template use File | Settings | File Templates.
 */
public interface StatisticalTest {
    public abstract Pair<Double,Double> getStatisticAndPValue();
    public abstract String getStatisticName();
}

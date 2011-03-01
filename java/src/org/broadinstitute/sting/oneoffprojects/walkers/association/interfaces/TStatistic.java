package org.broadinstitute.sting.oneoffprojects.walkers.association.interfaces;

import org.broadinstitute.sting.utils.collections.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 12:58:16 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TStatistic {
    public abstract double getTStatistic();
    public abstract double getDegreesOfFreedom();
}
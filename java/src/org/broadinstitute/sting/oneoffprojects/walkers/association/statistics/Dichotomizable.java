package org.broadinstitute.sting.oneoffprojects.walkers.association.statistics;

import org.broadinstitute.sting.utils.collections.Pair;

import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: 4/5/11
 * Time: 11:35 AM
 * To change this template use File | Settings | File Templates.
 */
public interface Dichotomizable {
    public abstract Pair<Collection<Number>,Collection<Number>> getDichotomizedData();
}

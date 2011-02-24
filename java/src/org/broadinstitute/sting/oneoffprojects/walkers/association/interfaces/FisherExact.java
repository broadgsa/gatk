package org.broadinstitute.sting.oneoffprojects.walkers.association.interfaces;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 12:58:56 PM
 * To change this template use File | Settings | File Templates.
 */
public interface FisherExact {
    public abstract int[][] getCounts();
    public abstract String[] getRowNames();
    public abstract String[] getColNames();
}
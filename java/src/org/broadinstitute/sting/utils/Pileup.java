package org.broadinstitute.sting.utils;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 9:33:00 AM
 * To change this template use File | Settings | File Templates.
 */
public interface Pileup {
    GenomeLoc getLocation();
    char getRef();
    String getBases();
    String getQuals();
    String getPileupString();
    int size();
}

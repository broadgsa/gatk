package org.broadinstitute.sting.utils;

import java.util.ArrayList;

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

    // byte array
    byte[] getBases();
    byte[] getQuals();

    ArrayList<Byte> getBasesAsArrayList();
    ArrayList<Byte> getQualsAsArrayList();

    String getBasesAsString();
    String getQualsAsString();
    String getPileupString();
    int size();
}

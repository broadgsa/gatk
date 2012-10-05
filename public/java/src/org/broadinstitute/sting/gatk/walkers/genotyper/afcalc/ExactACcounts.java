package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import java.util.Arrays;

/**
* Created with IntelliJ IDEA.
* User: depristo
* Date: 10/5/12
* Time: 2:54 PM
* To change this template use File | Settings | File Templates.
*/ // a wrapper around the int array so that we can make it hashable
public final class ExactACcounts {
    private final int[] counts;
    private int hashcode = -1;

    public ExactACcounts(final int[] counts) {
        this.counts = counts;
    }

    public int[] getCounts() {
        return counts;
    }

    @Override
    public boolean equals(Object obj) {
        return (obj instanceof ExactACcounts) && Arrays.equals(getCounts(), ((ExactACcounts) obj).getCounts());
    }

    @Override
    public int hashCode() {
        if ( hashcode == -1 )
            hashcode = Arrays.hashCode(getCounts());
        return hashcode;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getCounts()[0]);
        for ( int i = 1; i < getCounts().length; i++ ) {
            sb.append("/");
            sb.append(getCounts()[i]);
        }
        return sb.toString();
    }
}

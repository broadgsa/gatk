package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:49:47 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReferenceOrderedDatum implements Comparable<ReferenceOrderedDatum> {
    protected String name;

    public ReferenceOrderedDatum(String name) {
        this.name = name;
    }

    public String getName() { return this.name; }

    public abstract void parseLine(final String[] parts);

    public abstract String toString();
    public abstract String toSimpleString();
    public abstract String repl();

    public abstract GenomeLoc getLocation();
    public int compareTo( ReferenceOrderedDatum that ) {
        return getLocation().compareTo(that.getLocation());
    }

    public final String getContig() { return getLocation().getContig(); }
    public final long getStart() { return getLocation().getStart(); }
    public final long getStop() { return getLocation().getStop(); }
}

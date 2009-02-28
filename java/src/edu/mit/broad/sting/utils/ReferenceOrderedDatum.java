package edu.mit.broad.sting.utils;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:49:47 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReferenceOrderedDatum {
    public ReferenceOrderedDatum() { }

    public abstract void parseLine(final String[] parts);

    public abstract String toString();
    public abstract String toSimpleString();

    public abstract String getContig();
    public abstract long getStart();
    public abstract long getStop();
}

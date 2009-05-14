package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:49:47 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class BasicReferenceOrderedDatum implements ReferenceOrderedDatum {
    protected String name;

    public BasicReferenceOrderedDatum(String name) {
        this.name = name;
    }

    public String getName() { return this.name; }

    public abstract boolean parseLine(final Object header, final String[] parts) throws IOException;

    public abstract String toString();
    public String toSimpleString() { return toString(); }
    public String repl() { return this.toString(); }

    public String delimiterRegex() {
        return "\t";
    }

    public abstract GenomeLoc getLocation();
    
    public int compareTo( ReferenceOrderedDatum that ) {
        return getLocation().compareTo(that.getLocation());
    }

    public Object initialize(final File source) throws FileNotFoundException {
        //System.out.printf("Initialize called with %s%n", source);
        return null;
    }
}
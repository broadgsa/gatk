package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;

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
public interface ReferenceOrderedDatum extends Comparable<ReferenceOrderedDatum>, HasGenomeLocation {
    public String getName();
    public boolean parseLine(final Object header, final String[] parts) throws IOException;
    public String toString();
    public String toSimpleString();
    public String repl();

    /**
     * Used by the ROD system to determine how to split input lines
     * @return Regex string delimiter separating fields
     */
    public String delimiterRegex();

    public GenomeLoc getLocation();
    public int compareTo( ReferenceOrderedDatum that );

    /**
     * Backdoor hook to read header, meta-data, etc. associated with the file.  Will be
     * called by the ROD system before streaming starts
     *
     * @param source source data file on disk from which this rod stream will be pulled
     * @return a header object that will be passed to parseLine command
     */
    public Object initialize(final File source) throws FileNotFoundException;
}

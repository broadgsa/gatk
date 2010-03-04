/*
 * Copyright (c) 2010.  The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks;

import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;

import java.io.File;
import java.util.Iterator;


/**
 * @author aaron
 *         <p/>
 *         Class RMDTrack
 *         <p/>
 *         the basics of what a reference metadata track must contain.
 */
public abstract class RMDTrack {

    // the basics of a track:
    private final Class type;           // our type
    private final String name;          // the name
    private final File file;            // the associated file we create the reader from

    /**
     * Create a track
     *
     * @param type the type of track, used for track lookup
     * @param name the name of this specific track
     * @param file the associated file, for reference or recreating the reader
     */
    protected RMDTrack(Class type, String name, File file) {
        this.type = type;
        this.name = name;
        this.file = file;
    }

    public Class getType() {
        return type;
    }

    public String getName() {
        return name;
    }

    public File getFile() {
        return file;
    }

    /**
     * @return how to get an iterator of the underlying data.  This is all a track has to support,
     *         but other more advanced tracks support the query interface
     */
    public abstract Iterator<GATKFeature> getIterator();

    /**
     * helper function for determining if we are the same track
     *
     * @param name the name to match
     * @param type the type to match
     *
     * @return true on a match, false if the name or type is different
     */
    public boolean matches(String name, String type) {
        return (name.equals(this.name) && type.equals(this.type.getSimpleName()));
    }

    /**
     * do we support the query interface?
     * @return true if we can be cast to the QueryableTrack interface 
     */
    public abstract boolean supportsQuery();
}

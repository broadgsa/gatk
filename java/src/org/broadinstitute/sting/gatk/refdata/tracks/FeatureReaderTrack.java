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

import org.broad.tribble.FeatureReader;
import org.broadinstitute.sting.gatk.refdata.utils.FeatureToGATKFeatureIterator;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;


/**
 * 
 * @author aaron 
 * 
 * Class FeatureReaderTrack
 *
 * A feature reader track, implementing the RMDTrack for tracks that are generated out of Tribble
 */
public class FeatureReaderTrack extends RMDTrack implements QueryableTrack {
    // our feature reader - allows queries
    private FeatureReader reader;

    /**
     * Create a track
     *
     * @param type the type of track, used for track lookup
     * @param name the name of this specific track
     * @param file the associated file, for reference or recreating the reader
     * @param reader the feature reader to use as the underlying data source
     */
    public FeatureReaderTrack(Class type, String name, File file, FeatureReader reader) {
        super(type, name, file);
        this.reader = reader;
    }

    /**
     * @return how to get an iterator of the underlying data.  This is all a track has to support,
     *         but other more advanced tracks support the query interface
     */
    @Override
    public Iterator<GATKFeature> getIterator() {
        try {
            return new FeatureToGATKFeatureIterator(reader.iterator(),this.getName());
        } catch (IOException e) {
            throw new StingException("Unable to read from file",e);
        }
    }

    /**
     * do we support the query interface?
     *
     * @return true
     */
    @Override
    public boolean supportsQuery() {
        return true;
    }

    @Override
    public Iterator<GATKFeature> query(GenomeLoc interval) throws IOException {
        return new FeatureToGATKFeatureIterator(reader.query(interval.getContig(),(int)interval.getStart(),(int)interval.getStop()),this.getName());
    }

    @Override
    public Iterator<GATKFeature> query(GenomeLoc interval, boolean contained) throws IOException {
        return new FeatureToGATKFeatureIterator(reader.query(interval.getContig(),(int)interval.getStart(),(int)interval.getStop(), contained),this.getName());
    }

    @Override
    public Iterator<GATKFeature> query(String contig, int start, int stop) throws IOException {
        return new FeatureToGATKFeatureIterator(reader.query(contig,start,stop),this.getName());
    }

    @Override
    public Iterator<GATKFeature> query(String contig, int start, int stop, boolean contained) throws IOException {
        return new FeatureToGATKFeatureIterator(reader.query(contig,start,stop, contained),this.getName());
    }

    @Override
    public void close() {
        try {
            reader.close();
        } catch (IOException e) {
            throw new StingException("Unable to close reader " + reader.toString(),e);
        }
        reader = null;
    }
}

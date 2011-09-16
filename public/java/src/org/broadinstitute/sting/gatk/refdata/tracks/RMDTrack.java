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

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureSource;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broad.tribble.source.PerformanceLoggingFeatureSource;
import org.broadinstitute.sting.gatk.refdata.utils.FeatureToGATKFeatureIterator;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;


/**
 * @author aaron
 *         <p/>
 *         Class RMDTrack
 *         <p/>
 *         the basics of what a reference metadata track must contain.
 */
public class RMDTrack {
    private final static Logger logger = Logger.getLogger(RMDTrackBuilder.class);

    // the basics of a track:
    private final Class type;           // our type
    private final String name;          // the name
    private final File file;            // the associated file we create the reader from

    // our feature reader - allows queries
    private FeatureSource reader;

    // our sequence dictionary, which can be null
    private final SAMSequenceDictionary dictionary;

    /**
     * Parser to use when creating/parsing GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;

    // our codec type
    private final FeatureCodec codec;

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
     * Create a track
     *
     * @param type the type of track, used for track lookup
     * @param name the name of this specific track
     * @param file the associated file, for reference or recreating the reader
     * @param reader the feature reader to use as the underlying data source
     * @param dict the sam sequence dictionary
     * @param codec the feature codec we use to decode this type
     */
    public RMDTrack(Class type, String name, File file, FeatureSource reader, SAMSequenceDictionary dict, GenomeLocParser genomeLocParser, FeatureCodec codec) {
        this.type = type;
        this.name = name;
        this.file = file;
        this.reader = reader;
        this.dictionary = dict;
        this.genomeLocParser = genomeLocParser;
        this.codec = codec;
    }

    /**
     * @return how to get an iterator of the underlying data.  This is all a track has to support,
     *         but other more advanced tracks support the query interface
     */
    public CloseableIterator<GATKFeature> getIterator() {
        try {
            return new FeatureToGATKFeatureIterator(genomeLocParser,reader.iterator(),this.getName());
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(getFile(), "Unable to read from file", e);
        }
    }

    public CloseableIterator<GATKFeature> query(GenomeLoc interval) throws IOException {
        CloseableTribbleIterator<Feature> iter = reader.query(interval.getContig(),interval.getStart(),interval.getStop());
        if ( RMDTrackBuilder.MEASURE_TRIBBLE_QUERY_PERFORMANCE )
            logger.warn("Query " + getName() + ":" + ((PerformanceLoggingFeatureSource)reader).getPerformanceLog());
        return new FeatureToGATKFeatureIterator(genomeLocParser, iter, this.getName());
    }

    public void close() {
        try {
            reader.close();
        } catch (IOException e) {
            throw new UserException.MalformedFile("Unable to close reader " + reader.toString(),e);
        }
        reader = null;
    }

    public FeatureSource getReader() {
        return reader;
    }

    /**
     * get the sequence dictionary from the track, if available
     * @return a SAMSequenceDictionary if available, null if unavailable
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return dictionary;
    }

    public Object getHeader() {
        return reader.getHeader();
    }

    public FeatureCodec getCodec() {
        return codec;
    }
}
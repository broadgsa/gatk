/*
 * Copyright (c) 2010 The Broad Institute
 *
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks.builders;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.index.LinearIndex;
import org.broad.tribble.index.LinearIndexCreator;
import org.broadinstitute.sting.gatk.refdata.tracks.FeatureReaderTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class TribbleRMDTrackBuilder
 *
 * This class keeps track of the available codecs, and knows how to put together a track of
 * that gets iterators from the FeatureReader using Tribble.
 */
public class TribbleRMDTrackBuilder extends PluginManager<FeatureCodec> implements RMDTrackBuilder {

    // the linear index extension
    private static final String linearIndexExtension = ".idx";

    /** Create a new plugin manager. */
    public TribbleRMDTrackBuilder() {
        super(FeatureCodec.class, "Codecs", "Codec");
    }

    /** @return a list of all available tracks we currently have access to create */
    @Override
    public Map<String, Class> getAvailableTrackNamesAndTypes() {
        Map<String, Class> classes = new HashMap<String, Class>();
        //for (String c : this.pluginsByName.keySet()) // TODO: Aaron uncomment these two lines when Tribble is live
        //     if (!c.contains("SNP")) classes.put(c,this.pluginsByName.get(c));
        return classes;
    }

    /**
     * create a RMDTrack of the specified type
     *
     * @param targetClass the target class of track
     * @param name        what to call the track
     * @param inputFile   the input file
     *
     * @return an instance of the track
     * @throws RMDTrackCreationException
     *          if we don't know of the target class or we couldn't create it
     */
    @Override
    public RMDTrack createInstanceOfTrack(Class targetClass, String name, File inputFile) throws RMDTrackCreationException {
        // make a feature reader
        FeatureReader reader;
        try {
            FeatureCodec codec = this.createByType(targetClass);

            // check to see if the input file has an index
            if (!(new File(inputFile.getAbsolutePath() + linearIndexExtension).canRead())) {
                LinearIndex index = createIndex(inputFile, codec);
                reader = new FeatureReader(inputFile,index, codec);
            }
            else {
                reader = new FeatureReader(inputFile,codec);
            }
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to create reader with file " + inputFile, e);
        } catch (IOException e) {
            throw new StingException("Unable to make the index file for " + inputFile, e);
        }
        // return a feature reader track
        return new FeatureReaderTrack(targetClass, name, inputFile, reader);
    }

    /**
     * create an index for the input file
     * @param inputFile the input file
     * @param codec the codec to use
     * @return a linear index for the specified type
     * @throws IOException if we cannot write the index file
     */
    private LinearIndex createIndex(File inputFile, FeatureCodec codec) throws IOException {
        LinearIndexCreator create = new LinearIndexCreator(inputFile, codec);
        return create.createIndex();
    }
}

/**
 * a fake Tribble track, used to test out the Tribble interface and feature codec detection
 */
class FakeTribbleTrack implements FeatureCodec {

    @Override
    public Feature decode(String s) {
        return null;
    }

    @Override
    public int headerLineCount(File file) {
        return 0;
    }
}
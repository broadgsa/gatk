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

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.index.linear.LinearIndexCreator;
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
 *
 * Here's an example run command to find SNPs 200 base pairs up and downstream of the target file.
 *
 * java -jar dist/GenomeAnalysisTK.jar \
 * -R /broad/1KG/reference/human_b36_both.fasta \
 * -L 1:1863 \
 * -L MT:16520 \
 * -db /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_129_b36.rod \
 * -dbw 200 \
 * -l INFO \
 * -T DbSNPWindowCounter
 */
public class TribbleRMDTrackBuilder extends PluginManager<FeatureCodec> implements RMDTrackBuilder {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(TribbleRMDTrackBuilder.class);


    // the linear index extension
    public static final String linearIndexExtension = ".idx";

    /** Create a new plugin manager. */
    public TribbleRMDTrackBuilder() {
        super(FeatureCodec.class, "Codecs", "Codec");
    }

    /** @return a list of all available tracks we currently have access to create */
    @Override
    public Map<String, Class> getAvailableTrackNamesAndTypes() {
        Map<String, Class> classes = new HashMap<String, Class>();
        for (String c : this.pluginsByName.keySet())
             classes.put(c,this.pluginsByName.get(c));
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
        reader = createFeatureReader(targetClass, inputFile);
        // return a feature reader track
        return new FeatureReaderTrack(targetClass, name, inputFile, reader);
    }

    /**
     * create a feature reader of the specified type
     * @param targetClass the target codec type
     * @param inputFile the input file to create the track from (of the codec type)
     * @return the FeatureReader instance
     */
    public FeatureReader createFeatureReader(Class targetClass, File inputFile) {
        FeatureReader reader = null;
        try {
            // check to see if the input file has an index
            if (requireIndex(inputFile)) {
                logger.warn("Creating Tribble Index for file " + inputFile);
                LinearIndex index = createIndex(inputFile, this.createByType(targetClass));
                reader = new FeatureReader(inputFile,index, this.createByType(targetClass));
            }
            else {
                reader = new FeatureReader(inputFile,this.createByType(targetClass));
            }
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to create reader with file " + inputFile, e);
        } catch (IOException e) {
            throw new StingException("Unable to make the index file for " + inputFile, e);
        }
        return reader;
    }

    /**
     * create an index for the input file
     * @param inputFile the input file
     * @param codec the codec to use
     * @return a linear index for the specified type
     * @throws IOException if we cannot write the index file
     */
    public static LinearIndex createIndex(File inputFile, FeatureCodec codec) throws IOException {
        LinearIndexCreator create = new LinearIndexCreator(inputFile, codec);
        
        // if we can write the index, we should, but if not just create it in memory
        File indexFile = new File(inputFile.getAbsoluteFile() + linearIndexExtension);
        if (indexFile.getParentFile().canWrite() && (!indexFile.exists() || indexFile.canWrite()))
            return create.createIndex();
        else {
            logger.info("Unable to write to location " + indexFile + " for index file, creating index in memory only");
            return create.createIndex(null);
        }

    }

    /**
     * this function checks if we need to make an index file. There are three cases:
     * 1. The index file doesn't exist; return true
     * 2. The index does exist, but is older than the file.  We delete the index and return true
     * 3. else return false;
     * @param inputFile the target file to make an index for
     * @return true if we need to create an index, false otherwise
     */
    public static boolean requireIndex(File inputFile) {
        // can we read the index? if not, create an index
        File indexFile = new File(inputFile.getAbsolutePath() + linearIndexExtension);
        if (!(indexFile.canRead())) return true;
        if (inputFile.lastModified() > indexFile.lastModified()) {
            logger.warn("Removing out of date (index file date older than target file ) index file " + indexFile);
            indexFile.delete();
            return true;
        }
        return false;
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
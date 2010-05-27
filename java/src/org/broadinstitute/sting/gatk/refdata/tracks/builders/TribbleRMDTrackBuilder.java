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
import org.broad.tribble.*;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.index.linear.LinearIndexCreator;
import org.broad.tribble.readers.BasicFeatureReader;
import org.broad.tribble.util.LineReader;
import org.broadinstitute.sting.gatk.refdata.tracks.FeatureReaderTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException;
import org.broadinstitute.sting.utils.file.FSLock;
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
            classes.put(c, this.pluginsByName.get(c));
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
        // return a feature reader track
        return new FeatureReaderTrack(targetClass, name, inputFile, createFeatureReader(targetClass, inputFile));
    }

    /**
     * create a feature reader of the specified type
     * @param targetClass the target codec type
     * @param inputFile the input file to create the track from (of the codec type)
     * @return the FeatureReader instance
     */
    public FeatureReader createFeatureReader(Class targetClass, File inputFile) {
        FeatureReader reader = null;
        if (inputFile.getAbsolutePath().endsWith(".gz"))
            reader = createBasicFeatureReaderNoAssumedIndex(targetClass, inputFile);
        else
            reader = getLinearFeatureReader(targetClass, inputFile);
        return reader;
    }

    /**
     * create a feature reader, without assuming there exists an index.  This code assumes the feature
     * reader of the appropriate type will figure out what the right index type is, and determine if it
     * exists.
     *
     * @param targetClass the codec class type
     * @param inputFile the file to load
     * @return a feature reader implementation
     */
    private BasicFeatureReader createBasicFeatureReaderNoAssumedIndex(Class targetClass, File inputFile) {
        // we might not know the index type, try loading with the default reader constructor
        logger.debug("Attempting to blindly load " + inputFile);
        try {
            return new BasicFeatureReader(inputFile.getAbsolutePath(),this.createByType(targetClass));
        } catch (IOException e) {
            throw new StingException("Unable to create feature reader from file " + inputFile);
        }
    }

    /**
     * create a linear feature reader, where we create the index ahead of time
     * @param targetClass the target class
     * @param inputFile the tribble file to parse
     * @return the input file as a FeatureReader
     */
    private FeatureReader getLinearFeatureReader(Class targetClass, File inputFile) {
        FeatureReader reader;
        try {
            Index index = loadIndex(inputFile, this.createByType(targetClass), true);
            reader = new BasicFeatureReader(inputFile.getAbsolutePath(), index, this.createByType(targetClass));
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
     * @param onDisk write the index to disk?
     * @return a linear index for the specified type
     * @throws IOException if we cannot write the index file
     */
    public synchronized static Index loadIndex(File inputFile, FeatureCodec codec, boolean onDisk) throws IOException {

        // our return index
        LinearIndex returnIndex = null;

        // create the index file name, locking on the index file name
        File indexFile = new File(inputFile.getAbsoluteFile() + linearIndexExtension);
        FSLock lock = new FSLock(indexFile);

        // acquire a lock on the file
        boolean obtainedLock = lock.lock();
        try {
            // if the file exists, and we can read it, load the index from disk
            if (indexFile.exists() && indexFile.canRead() && obtainedLock) {

                // check to see if the index file is out of date
                if (indexFile.lastModified() < inputFile.lastModified())
                    logger.warn("Tribble index file " + indexFile + " is older than the track file " + inputFile + ", this can lead to unexpected behavior");
                
                logger.info("Loading Tribble index from disk for file " + inputFile);
                return LinearIndex.createIndex(indexFile);
            }
            // else we need to create the index, and write it to disk if we can
            else
                return writeIndexToDisk(inputFile, codec, onDisk, indexFile, obtainedLock);            
        }
        finally {
            lock.unlock();
        }

    }

    /**
     * attempt to create the index, and to disk
     * @param inputFile the input file
     * @param codec the codec to use
     * @param onDisk if they asked for disk storage or now
     * @param indexFile the index file location
     * @param obtainedLock did we obtain the lock on the file?
     * @return the index object
     * @throws IOException
     */
    private static LinearIndex writeIndexToDisk(File inputFile, FeatureCodec codec, boolean onDisk, File indexFile, boolean obtainedLock) throws IOException {
        LinearIndexCreator create = new LinearIndexCreator(inputFile, codec);
        LinearIndex index = create.createIndex(null); // we don't want to write initially, so we pass in null

        // if the index doesn't exist, and we can write to the directory, and we got a lock: write to the disk
        if (indexFile.getParentFile().canWrite() &&
                (!indexFile.exists() || indexFile.canWrite()) &&
                onDisk &&
                obtainedLock) {
            logger.info("Creating Tribble Index on disk for file " + inputFile);
            index.write(indexFile);
            return index;
        }
        // we can't write it to disk, just store it in memory
        else {
            // if they wanted to write, let them know we couldn't
            if (onDisk) logger.info("Unable to write to " + indexFile + " for the index file, creating index in memory only");
            return index;
        }
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
    public int readHeader(LineReader reader) {
        return 0; // the basics
    }
}
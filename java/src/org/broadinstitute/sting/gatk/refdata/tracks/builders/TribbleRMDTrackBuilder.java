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

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broad.tribble.*;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexCreator;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.index.interval.IntervalIndexCreator;
import org.broad.tribble.index.linear.LinearIndexCreator;
import org.broad.tribble.source.BasicFeatureSource;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broadinstitute.sting.gatk.refdata.tracks.TribbleTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.file.FSLockWithShared;
import org.broadinstitute.sting.utils.file.FileSystemInabilityToLockException;

import java.io.*;
import java.util.*;


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
    public static final String indexExtension = ".idx";

    /** Create a new plugin manager. */
    public TribbleRMDTrackBuilder() {
        super(FeatureCodec.class, "Codecs", "Codec");
    }

    /** @return a list of all available track types we currently have access to create */
    @Override
    public Map<String, Class> getAvailableTrackNamesAndTypes() {
        return new HashMap<String, Class>(this.pluginsByName);
    }

    /** @return a list of all available track record types we currently have access to create */
    @Override
    public Map<String, Class> getAvailableTrackNamesAndRecordTypes() {
        Map<String, Class> classes = new HashMap<String, Class>();
        for (String name: this.pluginsByName.keySet()) {
            FeatureCodec codec = this.createByName(name);
            classes.put(name, codec.getFeatureType());
        }
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
        Pair<BasicFeatureSource, SAMSequenceDictionary> pair = createFeatureReader(targetClass, name, inputFile);
        if (pair == null) throw new StingException("Unable to make the feature reader for input file " + inputFile);
        return new TribbleTrack(targetClass, name, inputFile, pair.first, pair.second, createCodec(targetClass, name));
    }

    public Pair<BasicFeatureSource, SAMSequenceDictionary> createFeatureReader(Class targetClass, File inputFile) {
        return createFeatureReader(targetClass, "anonymous", inputFile);
    }

    /**
     * create a feature reader of the specified type
     * @param targetClass the target codec type
     * @param name the target name
     * @param inputFile the input file to create the track from (of the codec type)
     * @return the FeatureReader instance
     */
    public Pair<BasicFeatureSource, SAMSequenceDictionary> createFeatureReader(Class targetClass, String name, File inputFile) {
        Pair<BasicFeatureSource, SAMSequenceDictionary> pair;
        if (inputFile.getAbsolutePath().endsWith(".gz"))
            pair = createBasicFeatureSourceNoAssumedIndex(targetClass, name, inputFile);
        else
            pair = getLinearFeatureReader(targetClass, name, inputFile);
        return pair;
    }

    /**
     * create a feature reader, without assuming there exists an index.  This code assumes the feature
     * reader of the appropriate type will figure out what the right index type is, and determine if it
     * exists.
     *
     * @param targetClass the codec class type
     * @param name the name of the track
     * @param inputFile the file to load
     * @return a feature reader implementation
     */
    private Pair<BasicFeatureSource, SAMSequenceDictionary> createBasicFeatureSourceNoAssumedIndex(Class targetClass, String name, File inputFile) {
        // we might not know the index type, try loading with the default reader constructor
        logger.info("Attempting to blindly load " + inputFile + " as a tabix indexed file");
        try {
            return new Pair<BasicFeatureSource, SAMSequenceDictionary>(new BasicFeatureSource(inputFile.getAbsolutePath(), createCodec(targetClass, name)),null);
        } catch (IOException e) {
            throw new StingException("Unable to create feature reader from file " + inputFile);
        }
    }

    private FeatureCodec createCodec(Class targetClass, String name) {
        FeatureCodec codex = this.createByType(targetClass);
        if ( codex instanceof NameAwareCodec )
            ((NameAwareCodec)codex).setName(name);
        return codex;
    }

    /**
     * create a linear feature reader, where we create the index ahead of time
     * @param targetClass the target class
     * @param name the name of the codec
     * @param inputFile the tribble file to parse
     * @return the input file as a FeatureReader
     */
    private Pair<BasicFeatureSource, SAMSequenceDictionary> getLinearFeatureReader(Class targetClass, String name, File inputFile) {
        Pair<BasicFeatureSource, SAMSequenceDictionary> reader;
        try {
            Index index = loadIndex(inputFile, createCodec(targetClass, name), true);
            reader = new Pair<BasicFeatureSource, SAMSequenceDictionary>(new BasicFeatureSource(inputFile.getAbsolutePath(),
                                                                                                index,
                                                                                                createCodec(targetClass, name)),
                                                                                                sequenceSetToDictionary(index.getSequenceNames()));
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

        // create the index file name, locking on the index file name
        File indexFile = null;
        indexFile = new File(inputFile.getAbsoluteFile() + indexExtension);

        FSLockWithShared lock = new FSLockWithShared(indexFile);

        // acquire a lock on the file
        Index idx = null;
        if (indexFile.canRead())
            idx = attemptIndexFromDisk(inputFile, codec, indexFile, lock);

        // if we managed to make an index, return
        if (idx != null) return idx;

        // we couldn't read the file, or we fell out of the conditions above, continue on to making a new index
        return createNewIndex(inputFile, codec, onDisk, indexFile, lock);
    }

    /**
     * attempt to read the index from disk
     * @param inputFile the input file
     * @param codec the codec to read from
     * @param indexFile the index file itself
     * @param lock the lock file
     * @return an index, or null if we couldn't load one
     * @throws IOException if we fail for FS issues
     */
    protected static Index attemptIndexFromDisk(File inputFile, FeatureCodec codec, File indexFile, FSLockWithShared lock) throws IOException {
        boolean locked;
        try {
            locked = lock.sharedLock();
        }
        catch(FileSystemInabilityToLockException ex) {
            throw new StingException("Unexpected inability to lock exception", ex);
        }
        Index idx;
        try {
            if (!locked) // can't lock file
                idx = createIndexInMemory(inputFile, codec);
            else
                idx = loadFromDisk(inputFile, indexFile);
        } finally {
            if (locked) lock.unlock();
        }
        return idx;
    }

    /**
     * load the index from disk, checking for out of date indexes and old versions (both of which are deleted)
     * @param inputFile the input file
     * @param indexFile the input file, plus the index extension
     * @return an Index, or null if we're unable to load
     */
    public static Index loadFromDisk(File inputFile, File indexFile) {
        logger.info("Loading Tribble index from disk for file " + inputFile);
        Index index = IndexFactory.loadIndex(indexFile.getAbsolutePath());

        // check if the file is up-to date (filestamp and version check)
        if (index.isCurrentVersion() && indexFile.lastModified() > inputFile.lastModified())
            return index;
        else if (indexFile.lastModified() < inputFile.lastModified())
            logger.warn("Index file " + indexFile + " is out of date (index older than input file), deleting and updating the index file");
        else // we've loaded an old version of the index, we want to remove it <-- currently not used, but may re-enable
            logger.warn("Index file " + indexFile + " is out of date (old version), deleting and updating the index file");

        // however we got here, remove the index and return null
        boolean deleted = indexFile.delete();

        if (!deleted) logger.warn("Index file " + indexFile + " is out of date, but could not be removed; it will not be trusted (we'll try to rebuild an in-memory copy)");
        return null;
    }


    /**
     * attempt to create the index, and write it to disk
     * @param inputFile the input file
     * @param codec the codec to use
     * @param onDisk if they asked for disk storage or now
     * @param indexFile the index file location
     * @param lock the locking object
     * @return the index object
     * @throws IOException when unable to create the new index
     */
    private static Index createNewIndex(File inputFile, FeatureCodec codec, boolean onDisk, File indexFile, FSLockWithShared lock) throws IOException {
        Index index = createIndexInMemory(inputFile, codec);

        boolean locked = false; // could we exclusive lock the file?
        try {
            locked = lock.exclusiveLock();
            if (locked) {
                logger.info("Writing Tribble index to disk for file " + inputFile);
                LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
                index.write(stream);
                stream.close();
            }
            else // we can't write it to disk, just store it in memory, tell them this
                if (onDisk) logger.info("Unable to write to " + indexFile + " for the index file, creating index in memory only");
            return index;
        }
        catch(FileSystemInabilityToLockException ex) {
            throw new StingException("Unexpected inability to lock exception", ex);
        }
        finally {
            if (locked) lock.unlock();
        }

    }

    /**
     * create the index in memory, given the input file and feature codec
     * @param inputFile the input file
     * @param codec the codec
     * @return a LinearIndex, given the file location
     * @throws IOException when unable to create the index in memory
     */
    private static Index createIndexInMemory(File inputFile, FeatureCodec codec) throws IOException {
        // this can take a while, let them know what we're doing
        logger.info("Creating Tribble index in memory for file " + inputFile);
        IndexCreator creator;
        creator = new LinearIndexCreator(inputFile, codec, null);
        return creator.createIndex();
    }

    /**
     * convert a list of Strings into a sequence dictionary
     * @param contigList the contig list, in coordinate order, this is allowed to be null
     * @return a SAMSequenceDictionary, WITHOUT contig sizes
     */
    private static SAMSequenceDictionary sequenceSetToDictionary(LinkedHashSet<String> contigList) {
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        if (contigList == null) return dict;

        for (String name : contigList) {
            SAMSequenceRecord seq = new SAMSequenceRecord(name, 0);
            dict.addSequence(seq);
        }
        return dict;
    }

}

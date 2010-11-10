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
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.source.BasicFeatureSource;
import org.broad.tribble.source.CachingFeatureSource;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.AbstractGenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.file.FSLockWithShared;
import org.broadinstitute.sting.utils.file.FileSystemInabilityToLockException;
import org.broadinstitute.sting.utils.instrumentation.Sizeof;

import java.io.*;
import java.util.*;


/**
 * 
 * @author aaron 
 *                                           `
 * Class RMDTrackBuilder
 *
 * This class keeps track of the available codecs, and knows how to put together a track of
 * that gets iterators from the FeatureReader using Tribble.
 *
 */
public class RMDTrackBuilder extends PluginManager<FeatureCodec> {
    /**
     * our log, which we use to capture anything from this class
     */
    private final static Logger logger = Logger.getLogger(RMDTrackBuilder.class);

    // a constant we use for marking sequence dictionary entries in the Tribble index property list
    public static final String SequenceDictionaryPropertyPredicate = "DICT:";

    // the input strings we use to create RODs from
    private final List<RMDTriplet> inputs = new ArrayList<RMDTriplet>();

    private Map<String, Class> classes = null;

    // private sequence dictionary we use to set our tracks with
    private SAMSequenceDictionary dict = null;

    /**
     * Private genome loc parser to use when building out new locs.
     */
    private GenomeLocParser genomeLocParser;

    /** Create a new plugin manager. */
    public RMDTrackBuilder() {
        super(FeatureCodec.class, "Codecs", "Codec");
    }

    /**
     * Create a new RMDTrackBuilder, with dictionary and genomeLocParser predefined.
     * @param dict
     * @param genomeLocParser
     */
    public RMDTrackBuilder(SAMSequenceDictionary dict,GenomeLocParser genomeLocParser) {
        super(FeatureCodec.class, "Codecs", "Codec");
        setSequenceDictionary(dict,genomeLocParser);
    }

    /**
     *
     * @param dict the sequence dictionary to use as a reference for Tribble track contig length lookups
     */
    public void setSequenceDictionary(SAMSequenceDictionary dict,GenomeLocParser genomeLocParser) {
        this.dict = dict;
        this.genomeLocParser = genomeLocParser;
    }    

    /** @return a list of all available track types we currently have access to create */
    public Map<String, Class> getAvailableTrackNamesAndTypes() {
        classes = new HashMap<String, Class>();
        for (String name: this.pluginsByName.keySet()) {
            classes.put(name.toUpperCase(), pluginsByName.get(name));
        }
        return classes;
    }

    /** @return a list of all available track record types we currently have access to create */
    public Map<String, Class> getAvailableTrackNamesAndRecordTypes() {
        HashMap classToRecord = new HashMap<String, Class>();
        for (String name: this.pluginsByName.keySet()) {
            FeatureCodec codec = this.createByName(name);
            classToRecord.put(name, codec.getFeatureType());
        }
        return classToRecord;
    }

    /**
     * create a RMDTrack of the specified type
     *
     * @param genomeLocParser GenomeLocParser to use, if case track needs additional reference context.
     * @param targetClass the target class of track
     * @param name        what to call the track
     * @param inputFile   the input file
     *
     * @return an instance of the track
     * @throws RMDTrackCreationException
     *          if we don't know of the target class or we couldn't create it
     */
    public RMDTrack createInstanceOfTrack(Class targetClass, String name, File inputFile) throws RMDTrackCreationException {
        // return a feature reader track
        Pair<FeatureSource, SAMSequenceDictionary> pair = createFeatureReader(targetClass, name, inputFile);
        if (pair == null) throw new UserException.CouldNotReadInputFile(inputFile, "Unable to make the feature reader for input file");
        return new RMDTrack(targetClass, name, inputFile, pair.first, pair.second, genomeLocParser, createCodec(targetClass,name));
    }

    /**
     * create a tribble feature reader class, given the target class and the input file
     * @param targetClass the target class, of a Tribble Codec type
     * @param inputFile the input file, that corresponds to the feature type
     * @return a pair of <BasicFeatureSource, SAMSequenceDictionary>
     */
    public Pair<FeatureSource, SAMSequenceDictionary> createFeatureReader(Class targetClass, File inputFile) {
        return createFeatureReader(targetClass, "anonymous", inputFile);
    }

    /**
     * create a feature reader of the specified type
     * @param targetClass the target codec type
     * @param name the target name
     * @param inputFile the input file to create the track from (of the codec type)
     * @return the FeatureReader instance
     */
    public Pair<FeatureSource, SAMSequenceDictionary> createFeatureReader(Class targetClass, String name, File inputFile) {
        Pair<FeatureSource, SAMSequenceDictionary> pair;
        if (inputFile.getAbsolutePath().endsWith(".gz"))
            pair = createBasicFeatureSourceNoAssumedIndex(targetClass, name, inputFile);
        else
            pair = getFeatureSource(targetClass, name, inputFile);
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
    private Pair<FeatureSource, SAMSequenceDictionary> createBasicFeatureSourceNoAssumedIndex(Class targetClass, String name, File inputFile) {
        // we might not know the index type, try loading with the default reader constructor
        logger.info("Attempting to blindly load " + inputFile + " as a tabix indexed file");
        try {
            return new Pair<FeatureSource, SAMSequenceDictionary>(BasicFeatureSource.getFeatureSource(inputFile.getAbsolutePath(), createCodec(targetClass, name)),null);
        } catch (TribbleException e) {
            throw new UserException(e.getMessage(), e);
        }
    }

    /**
     * add a name to the codec, if it takes one
     * @param targetClass the class to create a codec for
     * @param name the name to assign this codec
     * @return the feature codec itself
     */
    private FeatureCodec createCodec(Class targetClass, String name) {
        FeatureCodec codex = this.createByType(targetClass);
        if ( codex instanceof NameAwareCodec )
            ((NameAwareCodec)codex).setName(name);
        if(codex instanceof ReferenceDependentFeatureCodec)
            ((ReferenceDependentFeatureCodec)codex).setGenomeLocParser(genomeLocParser);        
        return codex;
    }

    /**
     * create a feature source object given:
     * @param targetClass the target class
     * @param name the name of the codec
     * @param inputFile the tribble file to parse
     * @return the input file as a FeatureReader
     */
    private Pair<FeatureSource, SAMSequenceDictionary> getFeatureSource(Class targetClass, String name, File inputFile) {
        Pair<FeatureSource, SAMSequenceDictionary> reader;
        try {
            Index index = loadIndex(inputFile, createCodec(targetClass, name));
            try { logger.info(String.format("  Index for %s has size in bytes %d", inputFile, Sizeof.getObjectGraphSize(index))); }
            catch ( ReviewedStingException e) { }

            SAMSequenceDictionary dictFromIndex = getSequenceDictionaryFromProperties(index);

            // if we don't have a dictionary in the Tribble file, and we've set a dictionary for this builder, set it in the file if they match
            if (dictFromIndex.size() == 0 && dict != null) {
                File indexFile = Tribble.indexFile(inputFile);
                setIndexSequenceDictionary(index,dict,indexFile,true);
                dictFromIndex = getSequenceDictionaryFromProperties(index);
            }

            FeatureSource source = new BasicFeatureSource(inputFile.getAbsolutePath(), index, createCodec(targetClass, name));
            //source = new CachingFeatureSource(source, 100, 100000);
            reader = new Pair<FeatureSource, SAMSequenceDictionary>(source, dictFromIndex);
        } catch (TribbleException e) {
            throw new UserException(e.getMessage());
        } catch (IOException e) {
            throw new UserException("Unable to create the index file for " + inputFile, e);
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
    public synchronized Index loadIndex(File inputFile, FeatureCodec codec) throws IOException {
        // create the index file name, locking on the index file name
        File indexFile = Tribble.indexFile(inputFile);
        FSLockWithShared lock = new FSLockWithShared(indexFile);

        // acquire a lock on the file
        Index idx = null;
        if (indexFile.canRead())
            idx = attemptIndexFromDisk(inputFile, codec, indexFile, lock);

        // if we managed to make an index, return
        if (idx != null) return idx;

        // we couldn't read the file, or we fell out of the conditions above, continue on to making a new index
        return writeIndexToDisk(createIndexInMemory(inputFile, codec), indexFile, lock);
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
    protected Index attemptIndexFromDisk(File inputFile, FeatureCodec codec, File indexFile, FSLockWithShared lock) throws IOException {
        boolean locked;
        try {
            locked = lock.sharedLock();
        }
        catch(FileSystemInabilityToLockException ex) {
            throw new UserException.MissortedFile(inputFile, "Unexpected inability to lock exception", ex);
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
        if (index.isCurrentVersion() && indexFile.lastModified() >= inputFile.lastModified())
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
     * attempt to write the index to disk
     * @param index the index to write to disk
     * @param indexFile the index file location
     * @param lock the locking object
     * @return the index object
     * @throws IOException when unable to create the new index
     */
    private static Index writeIndexToDisk(Index index, File indexFile, FSLockWithShared lock) throws IOException {
        boolean locked = false; // could we exclusive lock the file?
        try {
            locked = lock.exclusiveLock(); // handle the case where we aren't locking anything
            if (locked) {
                logger.info("Writing Tribble index to disk for file " + indexFile);
                LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
                index.write(stream);
                stream.close();
            }
            else // we can't write it to disk, just store it in memory, tell them this
                logger.warn("Unable to write to " + indexFile + " for the index file, creating index in memory only");

            try { logger.info(String.format("  Index for %s has size in bytes %d", indexFile, Sizeof.getObjectGraphSize(index))); }
            catch ( ReviewedStingException e) { }

            return index;
        }
        catch(FileSystemInabilityToLockException ex) {
            throw new UserException.CouldNotCreateOutputFile(indexFile,"Unexpected inability to lock exception", ex);
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
    private Index createIndexInMemory(File inputFile, FeatureCodec codec) {
        // this can take a while, let them know what we're doing
        logger.info("Creating Tribble index in memory for file " + inputFile);
        Index idx = IndexFactory.createIndex(inputFile, codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
        setIndexSequenceDictionary(idx, dict, null, false);
        return idx;
    }

    /**
     * Returns a collection of track names that match the record type.
     * @param trackRecordType the record type specified in the @RMD annotation
     * @return a collection of available track record type names that match the record type
     */
    public Collection<String> getTrackRecordTypeNames(Class trackRecordType) {
        Set<String> names = new TreeSet<String>();
        if (trackRecordType == null)
            throw new IllegalArgumentException("trackRecordType value is null, please pass in an actual class object");

        for (Map.Entry<String, Class> availableTrackRecordType: getAvailableTrackNamesAndRecordTypes().entrySet()) {
            if (availableTrackRecordType.getValue() != null && trackRecordType.isAssignableFrom(availableTrackRecordType.getValue()))
                names.add(availableTrackRecordType.getKey());
        }
        return names;
    }

    /**
     * find the associated reference meta data
     *
     * @param argCollection the input arguments to the GATK.
     * @param engine the GATK engine to bind the tracks to
     *
     * @return a list of RMDTracks, one for each -B option
     */
    public List<RMDTrack> getReferenceMetaDataSources(AbstractGenomeAnalysisEngine engine, GATKArgumentCollection argCollection) {
        initializeConvenienceBindings(engine,argCollection);
        initializeFullBindings(engine,argCollection);
        // try and make the tracks given their requests
        return createRequestedTrackObjects();
    }

    private void initializeConvenienceBindings(AbstractGenomeAnalysisEngine engine, GATKArgumentCollection argCollection) {
        if (argCollection.DBSNPFile != null)
            inputs.add(new RMDTriplet(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, "dbsnp", argCollection.DBSNPFile));
    }

    /**
     * initialize our lists of bindings
     * @param engine The engine, used to populate tags.
     * @param argCollection input arguments to the GATK.
     */
    private void initializeFullBindings(AbstractGenomeAnalysisEngine engine,GATKArgumentCollection argCollection) {
        // NOTE: Method acts as a static.  Once the inputs have been passed once they are locked in.
        if (argCollection.RODBindings.size() == 0)
            return;

        for (String binding: argCollection.RODBindings) {
            if(engine != null) {
                if(engine.getTags(binding).size() != 2)
                    throw new UserException("Invalid syntax for -B (reference-ordered data) input flag.  " +
                                            "Please use the following syntax when providing reference-ordered " +
                                            "data: -B:<name>,<type> <filename>.");
                // Assume that if tags are present, those tags are name and type.
                // Name is always first, followed by type.
                List<String> parameters = engine.getTags(binding);
                String name = parameters.get(0);
                String type = parameters.get(1);
                inputs.add(new RMDTriplet(name,type,binding));
            }
        }
    }

    /**
     * create the requested track objects
     *
     * @return a list of the tracks, one for each of the requested input tracks
     */
    private List<RMDTrack> createRequestedTrackObjects() {
        // create of live instances of the tracks
        List<RMDTrack> tracks = new ArrayList<RMDTrack>();

        // create instances of each of the requested types
        for (RMDTriplet trip : inputs) {
            Class featureCodecClass = getAvailableTrackNamesAndTypes().get(trip.getType().toUpperCase());
            if (featureCodecClass == null)
                throw new UserException.BadArgumentValue("-B",trip.getType());
            tracks.add(createInstanceOfTrack(featureCodecClass, trip.getName(), new File(trip.getFile())));
        }
        return tracks;
    }


    // ---------------------------------------------------------------------------------------------------------
    // static functions to work with the sequence dictionaries of indexes
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * get the sequence dictionary from the track, if available.  If not, make it from the contig list that is always in the index
     * @param index the index file to use
     * @return a SAMSequenceDictionary if available, null if unavailable
     */
    public static SAMSequenceDictionary getSequenceDictionaryFromProperties(Index index) {
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (Map.Entry<String,String> entry : index.getProperties().entrySet()) {
            if (entry.getKey().startsWith(SequenceDictionaryPropertyPredicate))
                dict.addSequence(new SAMSequenceRecord(entry.getKey().substring(SequenceDictionaryPropertyPredicate.length() , entry.getKey().length()),
                                 Integer.valueOf(entry.getValue())));
        }
        return dict;
    }

    /**
     * create the sequence dictionary with the contig list; a backup approach
     * @param index the index file to use
     * @param dict the sequence dictionary to add contigs to
     * @return the filled-in sequence dictionary
     */
    private static SAMSequenceDictionary createSequenceDictionaryFromContigList(Index index, SAMSequenceDictionary dict) {
        LinkedHashSet<String> seqNames = index.getSequenceNames();
        if (seqNames == null) {
            return dict;
        }
        for (String name : seqNames) {
            SAMSequenceRecord seq = new SAMSequenceRecord(name, 0);
            dict.addSequence(seq);
        }
        return dict;
    }

    /**
     * set the sequence dictionary of the track.  This function checks that the contig listing of the underlying file is compatible.
     * (that each contig in the index is in the sequence dictionary).
     * @param dict the sequence dictionary
     * @param index the index file
     * @param indexFile the index file
     * @param rewriteIndex should we rewrite the index when we're done?
     *
     */
    public static void setIndexSequenceDictionary(Index index, SAMSequenceDictionary dict, File indexFile, boolean rewriteIndex) {
        if (dict == null) return;

        SAMSequenceDictionary currentDict = createSequenceDictionaryFromContigList(index, new SAMSequenceDictionary());
        // check that every contig in the RMD contig list is at least in the sequence dictionary we're being asked to set
        for (SAMSequenceRecord seq : currentDict.getSequences()) {
            if (dict.getSequence(seq.getSequenceName()) == null)
                throw new UserException.IncompatibleSequenceDictionaries("The sequence dictionary from the reference the GATK is running with is not compatible with the sequence " +
                                                      "dictionary in the Tribble file " + indexFile + ".  It doesn't contain the contig: " + seq.getSequenceName(),
                                                      "RMD Sequence Dictionary",
                                                      currentDict,
                                                      "Reference Sequence Dictionary",
                                                      dict);
            index.addProperty(SequenceDictionaryPropertyPredicate + dict.getSequence(seq.getSequenceName()).getSequenceName(), String.valueOf(dict.getSequence(seq.getSequenceName()).getSequenceLength()));
        }
        // re-write the index
        if (rewriteIndex) try {
            writeIndexToDisk(index,indexFile,new FSLockWithShared(indexFile));
        } catch (IOException e) {
            logger.warn("Unable to update index with the sequence dictionary for file " + indexFile + "; this will not effect your run of the GATK");
        }
    }
}

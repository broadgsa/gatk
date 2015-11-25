/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.refdata.tracks;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.log4j.Logger;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import org.broadinstitute.gatk.utils.SequenceDictionaryUtils;
import org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet.RMDStorageType;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.file.FSLockWithShared;
import org.broadinstitute.gatk.utils.instrumentation.Sizeof;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;


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
public class RMDTrackBuilder { // extends PluginManager<FeatureCodec> {
    /**
     * our log, which we use to capture anything from this class
     */
    private final static Logger logger = Logger.getLogger(RMDTrackBuilder.class);

    // private sequence dictionary we use to set our tracks with
    private final SAMSequenceDictionary dict;

    /**
     * Private genome loc parser to use when building out new locs.
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * Validation exclusions, for validating the sequence dictionary.
     */
    private ValidationExclusion.TYPE validationExclusionType;

    private final FeatureManager featureManager;

    // If true, do not attempt to create index files if they don't exist or are outdated, and don't
    // make any file lock acquisition calls on the index files.
    private final boolean disableAutoIndexCreation;

    // Map of file name -> new sample name used when performing on-the-fly sample renaming
    private final Map<String, String> sampleRenameMap;

    /**
     * Construct an RMDTrackerBuilder, allowing the user to define tracks to build after-the-fact.  This is generally
     * used when walkers want to directly manage the ROD system for whatever reason.  Before using this constructor,
     * please talk through your approach with the SE team.
     * @param dict Sequence dictionary to use.
     * @param genomeLocParser Location parser to use.
     * @param validationExclusionType Types of validations to exclude, for sequence dictionary verification.
     * @param disableAutoIndexCreation Do not auto-create index files, and do not use file locking when accessing index files.
     *                                 UNSAFE in general (because it causes us not to lock index files before reading them) --
     *                                 suitable only for test suite use.
     * @param sampleRenameMap Map of file name -> new sample name used when performing on-the-fly sample renaming
     */
    public RMDTrackBuilder(final SAMSequenceDictionary dict,
                           final GenomeLocParser genomeLocParser,
                           final ValidationExclusion.TYPE validationExclusionType,
                           final boolean disableAutoIndexCreation,
                           final Map<String, String> sampleRenameMap) {
        this.dict = dict;
        this.validationExclusionType = validationExclusionType;
        this.genomeLocParser = genomeLocParser;
        this.featureManager = new FeatureManager(ValidationExclusion.lenientVCFProcessing(validationExclusionType));
        this.disableAutoIndexCreation = disableAutoIndexCreation;
        this.sampleRenameMap = sampleRenameMap;
    }

    /**
     * Return the feature manager this RMDTrackBuilder is using the create tribble tracks
     *
     * @return
     */
    public FeatureManager getFeatureManager() {
        return featureManager;
    }

    /**
     * create a RMDTrack of the specified type
     *
     * @param fileDescriptor a description of the type of track to build.
     *
     * @return an instance of the track
     */
    public RMDTrack createInstanceOfTrack(final RMDTriplet fileDescriptor) {
        String name = fileDescriptor.getName();
        File inputFile = new File(fileDescriptor.getFile());

        FeatureManager.FeatureDescriptor descriptor = getFeatureManager().getByTriplet(fileDescriptor);
        if (descriptor == null)
            throw new UserException.BadArgumentValue("-B",fileDescriptor.getType());

        // return a feature reader track
        Pair<AbstractFeatureReader, SAMSequenceDictionary> pair;
        if (ArgumentTypeDescriptor.isCompressed(inputFile.toString()))
            pair = createTabixIndexedFeatureSource(descriptor, name, inputFile);
        else
            pair = getFeatureSource(descriptor, name, inputFile, fileDescriptor.getStorageType());
        if (pair == null) throw new UserException.CouldNotReadInputFile(inputFile, "Unable to make the feature reader for input file");

        validateVariantAgainstSequenceDictionary(name, descriptor.getName(), pair.first, pair.second);

        return new RMDTrack(descriptor.getCodecClass(), name, inputFile, pair.first, pair.second, genomeLocParser, createCodec(descriptor, name, inputFile));
    }

    /**
     * Validate the VCF dictionary against the sequence dictionary.
     *
     * @param name      the name of this specific track
     * @param descriptorName  the name of the feature
     * @param reader    the feature reader to use as the underlying data source
     * @param dict      the sam sequence dictionary
     */
    private void validateVariantAgainstSequenceDictionary(final String name, final String descriptorName, final AbstractFeatureReader reader, final SAMSequenceDictionary dict ) throws UserException {
        // only process if the variant is a VCF
        if ( name.equals("variant") && descriptorName.equals("VCF") ){
            if ( reader != null && dict != null && reader.getHeader() != null ){
                final List<VCFContigHeaderLine> contigs = ((VCFHeader) reader.getHeader()).getContigLines();
                if (contigs != null) {
                    // make the VCF dictionary from the contig header fields
                    final List<SAMSequenceRecord> vcfContigRecords = new ArrayList<SAMSequenceRecord>();
                    for (final VCFContigHeaderLine contig : contigs)
                        vcfContigRecords.add(contig.getSAMSequenceRecord());

                    // have VCF contig fields so can make a dictionary and compare it to the sequence dictionary
                    if (!vcfContigRecords.isEmpty()) {
                        final SAMSequenceDictionary vcfDictionary = new SAMSequenceDictionary(vcfContigRecords);
                        final SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary(dict.getSequences());

                        SequenceDictionaryUtils.validateDictionaries(logger, validationExclusionType, name, vcfDictionary, "sequence", sequenceDictionary, false, null);
                    }
                }
            }
        }
    }

    /**
     * Convenience method simplifying track creation.  Assume unnamed track based on a file rather than a stream.
     * @param codecClass Type of Tribble codec class to build.
     * @param inputFile Input file type to use.
     * @return An RMDTrack, suitable for accessing reference metadata.
     */
    public RMDTrack createInstanceOfTrack(Class codecClass, File inputFile) {
        final FeatureManager.FeatureDescriptor descriptor = getFeatureManager().getByCodec(codecClass);

        if (descriptor == null)
            throw new ReviewedGATKException("Unable to find type name for codec class " + codecClass.getName());

        return createInstanceOfTrack(new RMDTriplet("anonymous",descriptor.getName(),inputFile.getAbsolutePath(),RMDStorageType.FILE,new Tags()));
    }

    /**
     * create a feature reader, without assuming there exists an index.  This code assumes the feature
     * reader of the appropriate type will figure out what the right index type is, and determine if it
     * exists.
     *
     * @param descriptor the FeatureDescriptor describing the FeatureCodec we want to create
     * @param name the name of the track
     * @param inputFile the file to load
     * @return a feature reader implementation
     */
    private Pair<AbstractFeatureReader, SAMSequenceDictionary> createTabixIndexedFeatureSource(FeatureManager.FeatureDescriptor descriptor, String name, File inputFile) {
        // we might not know the index type, try loading with the default reader constructor
        logger.debug("Attempting to load " + inputFile + " as a tabix indexed file without validating it");
        try {
            // getFeatureReader will detect that it's Tabix
            return new Pair<>(AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), createCodec(descriptor, name, inputFile)), null);
        } catch (TribbleException e) {
            throw new UserException(e.getMessage(), e);
        }
    }

    /**
     * add a name to the codec, if it takes one
     * @param descriptor the class to create a codec for
     * @param name the name to assign this codec
     * @param inputFile input file that we will be decoding
     * @return the feature codec itself
     */
    private FeatureCodec createCodec(final FeatureManager.FeatureDescriptor descriptor, final String name, final File inputFile) {
        // The remappedSampleName will be null if either no on-the-fly sample renaming was requested,
        // or the user's sample rename map file didn't contain an entry for this file:
        final String remappedSampleName = sampleRenameMap != null ? sampleRenameMap.get(inputFile.getAbsolutePath()) : null;

        return featureManager.createCodec(descriptor, name, genomeLocParser, remappedSampleName);
    }

    /**
     * create a feature source object given:
     * @param descriptor the FeatureDescriptor describing the FeatureCodec we want to create
     * @param name the name of the codec
     * @param inputFile the tribble file to parse
     * @param storageType How the RMD is streamed into the input file.
     * @return the input file as a FeatureReader
     */
    private Pair<AbstractFeatureReader, SAMSequenceDictionary> getFeatureSource(FeatureManager.FeatureDescriptor descriptor,
                                                                        String name,
                                                                        File inputFile,
                                                                        RMDStorageType storageType) {
        // Feature source and sequence dictionary to use as the ultimate reference
        AbstractFeatureReader featureSource = null;
        SAMSequenceDictionary sequenceDictionary = null;

        // Detect whether or not this source should be indexed.
        boolean canBeIndexed = (storageType == RMDStorageType.FILE);

        if(canBeIndexed) {
            try {
                Index index = loadIndex(inputFile, createCodec(descriptor, name, inputFile));
                try { logger.info(String.format("  Index for %s has size in bytes %d", inputFile, Sizeof.getObjectGraphSize(index))); }
                catch (ReviewedGATKException e) { }

                sequenceDictionary = IndexDictionaryUtils.getSequenceDictionaryFromProperties(index);

                // if we don't have a dictionary in the Tribble file, and we've set a dictionary for this builder, set it in the file if they match
                if (sequenceDictionary.isEmpty() && dict != null) {
                    validateAndUpdateIndexSequenceDictionary(inputFile, index, dict);

                    if ( ! disableAutoIndexCreation ) {
                        File indexFile = Tribble.indexFile(inputFile);
                        try { // re-write the index
                            writeIndexToDisk(index,indexFile,new FSLockWithShared(indexFile));
                        } catch (IOException e) {
                            logger.warn("Unable to update index with the sequence dictionary for file " + indexFile + "; this will not affect your run of the GATK");
                        }
                    }

                    sequenceDictionary = IndexDictionaryUtils.getSequenceDictionaryFromProperties(index);
                }

                featureSource = AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), createCodec(descriptor, name, inputFile), index);
            }
            catch (TribbleException e) {
                throw new UserException(e.getMessage());
            }
            catch (IOException e) {
                throw new UserException("I/O error loading or writing tribble index file for " + inputFile.getAbsolutePath(), e);
            }
        }
        else {
            featureSource = AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), createCodec(descriptor, name, inputFile), false);
        }

        return new Pair<AbstractFeatureReader,SAMSequenceDictionary>(featureSource,sequenceDictionary);
    }

    /**
     * create an index for the input file
     * @param inputFile the input file
     * @param codec the codec to use
     * @return a linear index for the specified type
     * @throws IOException if we cannot write the index file
     */
    public synchronized Index loadIndex( final File inputFile, final FeatureCodec codec) throws IOException {
        final File indexFile = Tribble.indexFile(inputFile);
        final FSLockWithShared lock = new FSLockWithShared(indexFile);
        Index idx = null;

        // If the index file exists and is readable, attempt to load it from disk. We'll get null back
        // if a problem was discovered with the index file when it was inspected, and we'll get an
        // in-memory index back in the case where the index file could not be locked.
        if (indexFile.canRead()) {
            idx = disableAutoIndexCreation ? loadFromDisk(inputFile, indexFile)  // load without locking if we're in disableAutoIndexCreation mode
                                           : attemptToLockAndLoadIndexFromDisk(inputFile, codec, indexFile, lock);
        }

        // If we have an index, it means we either loaded it from disk without issue or we created an in-memory
        // index due to not being able to acquire a lock.
        if (idx != null) return idx;

        // We couldn't read the file, or we discovered a problem with the index file, so continue on to making a new index
        idx = createIndexInMemory(inputFile, codec);
        if ( ! disableAutoIndexCreation ) {
            writeIndexToDisk(idx, indexFile, lock);
        }
        return idx;
    }

    /**
     * Attempt to acquire a shared lock and then load the index from disk. Returns an in-memory index if
     * a lock could not be obtained. Returns null if a problem was discovered with the index file when it
     * was examined (eg., it was out-of-date).
     *
     * @param inputFile the input file
     * @param codec the codec to read from
     * @param indexFile the index file itself
     * @param lock the lock file
     * @return an index, or null if we couldn't load one
     * @throws IOException if we fail for FS issues
     */
    protected Index attemptToLockAndLoadIndexFromDisk( final File inputFile, final FeatureCodec codec, final File indexFile, final FSLockWithShared lock ) throws IOException {
        boolean locked = false;
        Index idx = null;

        try {
            locked = lock.sharedLock();

            if ( ! locked ) { // can't lock file
                logger.info(String.format("Could not acquire a shared lock on index file %s, falling back to using an in-memory index for this GATK run.",
                                          indexFile.getAbsolutePath()));
                idx = createIndexInMemory(inputFile, codec);
            }
            else {
                idx = loadFromDisk(inputFile, indexFile);
            }
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
    protected Index loadFromDisk( final File inputFile, final File indexFile ) {
        logger.debug("Loading Tribble index from disk for file " + inputFile);
        Index index = IndexFactory.loadIndex(indexFile.getAbsolutePath());

        // check if the file is up-to date (filestamp and version check)
        if (index.isCurrentVersion() && indexFile.lastModified() >= inputFile.lastModified())
            return index;
        else if (indexFile.lastModified() < inputFile.lastModified())
            logger.warn("Index file " + indexFile + " is out of date (index older than input file), " +
                        (disableAutoIndexCreation ? "falling back to an in-memory index" : "deleting and updating the index file"));
        else // we've loaded an old version of the index, we want to remove it <-- currently not used, but may re-enable
            logger.warn("Index file " + indexFile + " is out of date (old version), " +
                        (disableAutoIndexCreation ? "falling back to an in-memory index" : "deleting and updating the index file"));

        if ( ! disableAutoIndexCreation ) {
            boolean deleted = indexFile.delete();
            if (!deleted) logger.warn("Index file " + indexFile + " is out of date, but could not be removed; it will not be trusted (we'll try to rebuild an in-memory copy)");
        }

        return null;
    }


    /**
     * attempt to write the index to disk
     * @param index the index to write to disk
     * @param indexFile the index file location
     * @param lock the locking object
     * @throws IOException when unable to create the new index
     */
    private void writeIndexToDisk( final Index index, final File indexFile, final FSLockWithShared lock ) throws IOException {
        if ( disableAutoIndexCreation ) {
            return;
        }

        boolean locked = false;

        try {
            locked = lock.exclusiveLock();

            if (locked) {
                logger.info("Writing Tribble index to disk for file " + indexFile);
                LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
                index.write(stream);
                stream.close();
            }
            else // we can't write it to disk, just store it in memory, tell them this
                logger.warn("Unable to write to " + indexFile + " for the index file, creating index in memory only");

            try { logger.info(String.format("  Index for %s has size in bytes %d", indexFile, Sizeof.getObjectGraphSize(index))); }
            catch ( ReviewedGATKException e) { }
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
    protected Index createIndexInMemory(File inputFile, FeatureCodec codec) {
        // this can take a while, let them know what we're doing
        logger.debug("Creating Tribble index in memory for file " + inputFile);
        Index idx = IndexFactory.createDynamicIndex(inputFile, codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
        validateAndUpdateIndexSequenceDictionary(inputFile, idx, dict);
        return idx;
    }

    /**
     * set the sequence dictionary of the track.  This function checks that the contig listing of the underlying file is compatible.
     * (that each contig in the index is in the sequence dictionary).
     * @param inputFile for proper error message formatting.
     * @param dict the sequence dictionary
     * @param index the index file
     */
    public void validateAndUpdateIndexSequenceDictionary(final File inputFile, final Index index, final SAMSequenceDictionary dict) {
        if (dict == null) throw new ReviewedGATKException("BUG: dict cannot be null");

        // check that every contig in the RMD contig list is at least in the sequence dictionary we're being asked to set
        final SAMSequenceDictionary currentDict = IndexDictionaryUtils.createSequenceDictionaryFromContigList(index, new SAMSequenceDictionary());
        validateTrackSequenceDictionary(inputFile.getAbsolutePath(), currentDict, dict);

        // actually update the dictionary in the index
        IndexDictionaryUtils.setIndexSequenceDictionary(index, dict);
    }

    public void validateTrackSequenceDictionary(final String trackName,
                                                final SAMSequenceDictionary trackDict,
                                                final SAMSequenceDictionary referenceDict ) {
        IndexDictionaryUtils.validateTrackSequenceDictionary(trackName, trackDict, referenceDict, validationExclusionType);
    }
}

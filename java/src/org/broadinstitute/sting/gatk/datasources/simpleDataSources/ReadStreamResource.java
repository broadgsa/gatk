/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.picard.sam.SamFileHeaderMerger;

import java.util.*;
import java.io.File;

/**
 * Represents a single stream of read data.  Used to represent the state of the stream and determine
 * whether the state of this resource is such that it can field the desired query.
 * @author hanna
 * @version 0.1
 */
class ReadStreamResource {
    final static boolean eagerDecode = true;

    /**
     * Do all the constituent components of this ReadStreamResource have indices?
     * In general, BAM files without indices are not supported, but in a few specific
     * cases we do allow this for the Picard pipeline.
     * @return true if all BAM files have indices; false otherwise.
     */
    protected boolean hasIndex() {
        for(SAMFileReader reader: readStreamPointer.getHeaderMerger().getReaders()) {
            if(!reader.hasIndex())
                return false;
        }
        return true;
    }

    protected final boolean hasIndex;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(ReadStreamPointer.class);

    /**
     * The (possibly merged) header for the input fileset.
     */
    private final SAMFileHeader header;

    /**
     * A pointer to the current location of the file.
     */
    private ReadStreamPointer readStreamPointer = null;

    /**
     * A mapping from original input file to merged read group record ids
     */
    private Map<File, Set<String>> fileToReadGroupIdMap = null;

    public ReadStreamResource( Reads sourceInfo ) {
        SamFileHeaderMerger headerMerger = createHeaderMerger(sourceInfo, SAMFileHeader.SortOrder.coordinate);

        this.header = headerMerger.getMergedHeader();

        boolean indexPresent = true;
        for(SAMFileReader reader: headerMerger.getReaders()) {
            if(!reader.hasIndex())
                indexPresent = false;
        }
        hasIndex = indexPresent;

        if(hasIndex)
            readStreamPointer = new MappedReadStreamPointer(sourceInfo, headerMerger);
        else
            readStreamPointer = new EntireReadStreamPointer(sourceInfo, headerMerger);
    }

    /**
     * Gets the header information for the read stream.
     * @return Header information for the read stream.
     */
    public SAMFileHeader getHeader() {
        return header;
    }
    
    /**
     * Returns Reads data structure containing information about the reads data sources as well as
     * information about how they are downsampled, sorted, and filtered
     * @return the Reads object
     */
    public Reads getReadsInfo() { return readStreamPointer.getReadsInfo(); }
    
    /** 
     * Returns header merger: a class that keeps the mapping between original read groups and read groups
     * of the merged stream; merger also provides access to the individual file readers (and hence headers
     * too) maintained by the system. 
     * @return the header merger
     */
   public SamFileHeaderMerger getHeaderMerger() { return readStreamPointer.getHeaderMerger(); }

    public boolean canAccessSegmentEfficiently(DataStreamSegment segment) {
        return readStreamPointer.canAccessSegmentEfficiently(segment);
    }

    public void close() {
        readStreamPointer.close();
    }

    public void destroy( StingSAMIterator iterator ) {
        readStreamPointer.destroy(iterator);
    }

    public StingSAMIterator getReadsContainedBy( DataStreamSegment segment ) {
        if( readStreamPointer instanceof MappedReadStreamPointer && segment instanceof UnmappedStreamSegment )
            readStreamPointer = ((MappedReadStreamPointer)readStreamPointer).toUnmappedReadStreamPointer();
        return readStreamPointer.getReadsContainedBy(segment);
    }

    public StingSAMIterator getReadsOverlapping( DataStreamSegment segment ) {
        return readStreamPointer.getReadsOverlapping(segment);
    }

    public Map<File, Set<String>> getFileToReadGroupIdMapping() {
        return fileToReadGroupIdMap;
    }

    /**
     * A private function that, given the internal file list, generates a merging construct for
     * all available files.
     * @param reads source information about the reads.
     * @param SORT_ORDER sort order for the reads.
     * @return a list of SAMFileReaders that represent the stored file names
     * @throws SimpleDataSourceLoadException if the file cannot be opened.
     */
    private SamFileHeaderMerger createHeaderMerger( Reads reads, SAMFileHeader.SortOrder SORT_ORDER )
            throws SimpleDataSourceLoadException {

        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        List<SAMFileReader> lst = new ArrayList<SAMFileReader>();
        Map<File, SAMFileReader> fileToReaderMap = new HashMap<File, SAMFileReader>();
        for (File f : reads.getReadsFiles()) {
            SAMFileReader reader = new SAMFileReader(f, eagerDecode);
            fileToReaderMap.put(f, reader);
            reader.setValidationStringency(reads.getValidationStringency());

            final SAMFileHeader header = reader.getFileHeader();
            logger.debug(String.format("Sort order is: " + header.getSortOrder()));

            if (reader.getFileHeader().getReadGroups().size() < 1) {
                //logger.warn("Setting header in reader " + f.getName());
                SAMReadGroupRecord rec = new SAMReadGroupRecord(f.getName());
                rec.setLibrary(f.getName());
                rec.setSample(f.getName());

                reader.getFileHeader().addReadGroup(rec);
            }

            lst.add(reader);
        }

        // create the header merger
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(lst,SORT_ORDER,true);

        // populate the file -> read group mapping
        fileToReadGroupIdMap = new HashMap<File, Set<String>>();
        for (Map.Entry<File, SAMFileReader> entry : fileToReaderMap.entrySet()) {

            Set<String> readGroups = new HashSet<String>(5);

            for (SAMReadGroupRecord g : entry.getValue().getFileHeader().getReadGroups()) {
                if (headerMerger.hasReadGroupCollisions()) {
                    // Check if there were read group clashes.
                    // If there were, use the SamFileHeaderMerger to translate from the
                    // original read group id to the read group id in the merged stream
                    readGroups.add(headerMerger.getReadGroupId(entry.getValue(), g.getReadGroupId()));
                } else {
                    // otherwise, pass through the unmapped read groups since this is what Picard does as well
                    readGroups.add(g.getReadGroupId());
                }
            }

            fileToReadGroupIdMap.put(entry.getKey(), readGroups);
        }

        return headerMerger;
    }
}

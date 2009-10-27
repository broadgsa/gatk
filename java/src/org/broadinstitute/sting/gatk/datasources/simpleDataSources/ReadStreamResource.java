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
import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.picard.sam.SamFileHeaderMerger;

import java.util.List;
import java.util.ArrayList;
import java.io.File;

/**
 * Represents a single stream of read data.  Used to represent the state of the stream and determine
 * whether the state of this resource is such that it can field the desired query.
 * @author hanna
 * @version 0.1
 */
class ReadStreamResource {
    final static boolean eagerDecode = true;

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

    public ReadStreamResource( Reads sourceInfo ) {
        SamFileHeaderMerger headerMerger = createHeaderMerger(sourceInfo, SAMFileHeader.SortOrder.coordinate);

        this.header = headerMerger.getMergedHeader();
        readStreamPointer = new MappedReadStreamPointer(sourceInfo, headerMerger);
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
     * @return
     */
    public Reads getReadsInfo() { return readStreamPointer.getReadsInfo(); }
    
    /** 
     * Returns header merger: a class that keeps the mapping between original read groups and read groups
     * of the merged stream; merger also provides access to the individual file readers (and hence headers
     * too) maintained by the system. 
     * @return
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


    public StingSAMIterator getReadsOverlapping( MappedStreamSegment segment ) {
        return readStreamPointer.getReadsOverlapping(segment);
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
        for (File f : reads.getReadsFiles()) {
            SAMFileReader reader = new SAMFileReader(f, eagerDecode);
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
        return new SamFileHeaderMerger(lst,SORT_ORDER,true);
    }
}

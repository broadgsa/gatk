/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.SamReader.Indexing;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * User: hanna
 * Date: Jun 11, 2009
 * Time: 9:35:31 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Pass specified reads into the given walker.
 */

public class ArtificialSAMFileReader implements SamReader, Indexing {

    /**
     * The reader of SamRecords
     */
    private SamReader reader;

    /**
     * The parser, for GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * Backing data store of reads.
     */
    private final List<SAMRecord> reads;

    /**
     * Input/custom SAM file header
     */
    private SAMFileHeader customHeader = null;

    /**
     * Construct an artificial SAM file reader.
     *
     * @param sequenceDictionary sequence dictionary used to initialize our GenomeLocParser
     * @param reads Reads to use as backing data source.
     */
    public ArtificialSAMFileReader(SAMSequenceDictionary sequenceDictionary,SAMRecord... reads) {
        final SamInputResource samInputResource = SamInputResource.of(createEmptyInputStream());
        reader = SamReaderFactory.makeDefault().open(samInputResource);
        this.genomeLocParser = new GenomeLocParser(sequenceDictionary);
        this.reads = Arrays.asList(reads);
    }

    /**
     * Construct an artificial SAM file reader with the given SAM file header
     *
     * @param customHeader Header that should be returned by calls to getFileHeader() on this reader
     * @param reads Reads to use as backing data source.
     */
    public ArtificialSAMFileReader( SAMFileHeader customHeader, SAMRecord... reads ) {
        final SamInputResource samInputResource = SamInputResource.of(createEmptyInputStream());
        reader = SamReaderFactory.makeDefault().open(samInputResource);

        this.customHeader = customHeader;
        this.genomeLocParser = new GenomeLocParser(customHeader.getSequenceDictionary());
        this.reads = Arrays.asList(reads);
    }

    @Override
    public String getResourceDescription() {
        return this.toString();
    }

    @Override
    public boolean hasIndex() {
        return this.reader.hasIndex();
    }

    @Override
    public Indexing indexing() {
        return this;
    }

    @Override
    public BrowseableBAMIndex getBrowseableIndex() {
        BAMIndex index = this.getIndex();
        if(!(index instanceof BrowseableBAMIndex)) {
            throw new SAMException("Cannot return index: index created by BAM is not browseable.");
        } else {
            return BrowseableBAMIndex.class.cast(index);
        }
    }

    @Override
    public boolean hasBrowseableIndex() {
        return this.hasIndex() && this.getIndex() instanceof BrowseableBAMIndex;
    }

    @Override
    public BAMIndex getIndex() {
        throw new UnsupportedOperationException();
    }

    @Override
    public SAMRecordIterator iterator() {
        return new SAMRecordIterator() {
            private final Iterator<SAMRecord> iterator = reads.iterator();
            public boolean hasNext() { return iterator.hasNext(); }
            public SAMRecord next() { return iterator.next(); }
            public void close() {}
            public void remove() { iterator.remove(); }
            public SAMRecordIterator assertSorted(SAMFileHeader.SortOrder sortOrder) { return this; }
        };
    }

    /**
     * Iterate through the the file.
     *
     * @param chunks List of chunks for which to retrieve data.
     * @return An iterator.
     */
    @Override
    public SAMRecordIterator iterator(SAMFileSpan chunks) {
        return new SamReader.AssertingIterator(this.reader.iterator());
    }

    public SAMRecordIterator query(final String sequence, final int start, final int end, final boolean contained) {
        GenomeLoc region = genomeLocParser.createGenomeLoc(sequence, start, end);
        List<SAMRecord> coveredSubset = new ArrayList<>();

        for( SAMRecord read: reads ) {
            GenomeLoc readPosition = genomeLocParser.createGenomeLoc(read);
            if( contained && region.containsP(readPosition) ) coveredSubset.add(read);
            else if( !contained && readPosition.overlapsP(region) ) coveredSubset.add(read);
        }

        final Iterator<SAMRecord> iterator = coveredSubset.iterator();
        return new SAMRecordIterator() {
            public boolean hasNext() { return iterator.hasNext(); }
            public SAMRecord next() { return iterator.next(); }
            public void close() {}
            public void remove() { iterator.remove(); }
            public SAMRecordIterator assertSorted(SAMFileHeader.SortOrder sortOrder) { return this; }
        };
    }

    @Override
    public SAMRecordIterator queryOverlapping(final String sequence, final int start, final int end) {
        return this.query(sequence, start, end, false);
    }

    @Override
    public SAMRecordIterator queryContained(final String sequence, final int start, final int end) {
        return this.query(sequence, start, end, true);
    }

    @Override
    public SAMRecordIterator query(final QueryInterval[] intervals, final boolean contained) {
        return new AssertingIterator(this.reader.query(intervals, contained));
    }

    @Override
    public SAMRecordIterator queryOverlapping(final QueryInterval[] intervals) {
        return this.query(intervals, false);
    }

    @Override
    public SAMRecordIterator queryContained(final QueryInterval[] intervals) {
        return this.query(intervals, true);
    }

    @Override
    public SAMRecordIterator queryUnmapped() {
        return new AssertingIterator(this.reader.queryUnmapped());
    }

    @Override
    public SAMRecordIterator queryAlignmentStart(final String sequence, final int start) {
        return new AssertingIterator(this.reader.queryAlignmentStart(sequence, start));
    }

    @Override
    public SAMRecord queryMate(final SAMRecord rec) {
        if(!rec.getReadPairedFlag()) {
            throw new IllegalArgumentException("queryMate called for unpaired read.");
        } else if(rec.getFirstOfPairFlag() == rec.getSecondOfPairFlag()) {
            throw new IllegalArgumentException("SAMRecord must be either first and second of pair, but not both.");
        } else {
            boolean firstOfPair = rec.getFirstOfPairFlag();
            SAMRecordIterator it;
            if(rec.getMateReferenceIndex() == -1) {
                it = this.queryUnmapped();
            } else {
                it = this.queryAlignmentStart(rec.getMateReferenceName(), rec.getMateAlignmentStart());
            }

            try {
                SAMRecord mateRec = null;

                while(true) {
                    SAMRecord next;
                    while(it.hasNext()) {
                        next = it.next();
                        if(!next.getReadPairedFlag()) {
                            if(rec.getReadName().equals(next.getReadName())) {
                                throw new SAMFormatException("Paired and unpaired reads with same name: " + rec.getReadName());
                            }
                        } else {
                            if(firstOfPair) {
                                if(next.getFirstOfPairFlag()) {
                                    continue;
                                }
                            } else if(next.getSecondOfPairFlag()) {
                                continue;
                            }

                            if(rec.getReadName().equals(next.getReadName())) {
                                if(mateRec != null) {
                                    throw new SAMFormatException("Multiple SAMRecord with read name " + rec.getReadName() + " for " + (firstOfPair?"second":"first") + " end.");
                                }

                                mateRec = next;
                            }
                        }
                    }

                    next = mateRec;
                    return next;
                }
            } finally {
                it.close();
            }
        }
    }

    @Override
    public SAMFileSpan getFilePointerSpanningReads() {
        return this.reader.indexing().getFilePointerSpanningReads();
    }

    @Override
    public void close() throws IOException{
        if(this.reader != null) {
            this.reader.close();
        }

        this.reader = null;
    }

    @Override
    public Type type() {
        return this.reader.type();
    }

    @Override
    public SAMFileHeader getFileHeader() {
        return customHeader != null ? customHeader : this.reader.getFileHeader();
    }

    /**
     * Builds an empty input stream for faking out the sam file reader.
     * Derive it from a string so that, in the future, it might be possible
     * to fake the text of a sam file from samtools output, et.c
     * @return Stream that returns no characters.
     */
    private static InputStream createEmptyInputStream() {
        try {
            byte[] byteArray = "".getBytes("ISO-8859-1");
            return new ByteArrayInputStream(byteArray);
        }
        catch( Exception ex ) {
            throw new ReviewedGATKException("Unable to build empty input stream",ex);
        }
    }
}

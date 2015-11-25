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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.ByteArrayInputStream;
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

public class ArtificialSAMFileReader extends SAMFileReader {
    /**
     * The parser, for GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * Backing data store of reads.
     */
    private final List<SAMRecord> reads;

    private SAMFileHeader customHeader = null;

    /**
     * Construct an artificial SAM file reader.
     * @param sequenceDictionary sequence dictionary used to initialize our GenomeLocParser
     * @param reads Reads to use as backing data source.
     */
    public ArtificialSAMFileReader(SAMSequenceDictionary sequenceDictionary,SAMRecord... reads) {
        super( createEmptyInputStream(),true );
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
        super(createEmptyInputStream(),true);

        this.customHeader = customHeader;
        this.genomeLocParser = new GenomeLocParser(customHeader.getSequenceDictionary());
        this.reads = Arrays.asList(reads);
    }


    @Override
    public SAMFileHeader getFileHeader() {
        if ( customHeader != null ) {
            return customHeader;
        }

        return super.getFileHeader();
    }

    /**
     * @{inheritDoc}
     */
    @Override
    public SAMRecordIterator query(final String sequence, final int start, final int end, final boolean contained) {
        GenomeLoc region = genomeLocParser.createGenomeLoc(sequence, start, end);
        List<SAMRecord> coveredSubset = new ArrayList<SAMRecord>();

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

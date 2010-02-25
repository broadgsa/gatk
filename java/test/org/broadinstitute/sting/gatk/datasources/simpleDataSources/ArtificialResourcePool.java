package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.sam.ArtificialSAMIterator;
import org.broadinstitute.sting.utils.sam.ArtificialSAMQueryIterator;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;

import java.util.Collections;
import java.io.File;


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

/**
 * use this to inject into SAMDataSource for testing
 */
public class ArtificialResourcePool extends SAMResourcePool {
    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.SILENT;

    // the header
    private SAMFileHeader header;
    private ArtificialSAMIterator iterator;

    /**
     * Track the iterator to see whether it's venturing into unmapped reads for the first
     * time.  If so, query straight there.  Only works for query iterators.
     *
     * TODO: Clean up.
     */
    private boolean intoUnmappedReads = false;

    public ArtificialResourcePool( SAMFileHeader header, ArtificialSAMIterator iterator ) {
        super( new Reads(Collections.<File>emptyList()) );
        this.header = header;
        this.iterator = iterator;
    }

    @Override
    public StingSAMIterator iterator( DataStreamSegment segment ) {
        if (segment instanceof MappedStreamSegment && iterator instanceof ArtificialSAMQueryIterator) {
            ArtificialSAMQueryIterator queryIterator = (ArtificialSAMQueryIterator)iterator;
            MappedStreamSegment mappedSegment = (MappedStreamSegment)segment;
            GenomeLoc bounds = mappedSegment.locus;
            if (!this.queryOverlapping) {
                queryIterator.queryContained(bounds.getContig(), (int)bounds.getStart(), (int)bounds.getStop());
            } else {
                queryIterator.queryOverlapping(bounds.getContig(), (int)bounds.getStart(), (int)bounds.getStop());
            }
            return queryIterator;
        }
        else if (segment instanceof UnmappedStreamSegment) {
            if( !intoUnmappedReads ) {
                if( iterator instanceof ArtificialSAMQueryIterator ) {
                    ArtificialSAMQueryIterator queryIterator = (ArtificialSAMQueryIterator)iterator;
                    queryIterator.queryUnmappedReads();
                }
                intoUnmappedReads = true;
            }
            return new BoundedReadIterator(iterator,((UnmappedStreamSegment)segment).size);
        }
        else
            throw new StingException("Unsupported segment type passed to test");
    }

    /**
     * get the merged header
     *
     * @return the merged header
     */
    public SAMFileHeader getHeader() {
       return this.header;
    }
}

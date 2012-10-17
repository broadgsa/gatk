/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Simple wrapper class that multiplexes multiple ArtificialSingleSampleReadStreams into a single stream of reads
 *
 * @author David Roazen
 */
public class ArtificialMultiSampleReadStream implements Iterable<SAMRecord> {

    private Collection<ArtificialSingleSampleReadStream> perSampleArtificialReadStreams;
    private MergingSamRecordIterator mergingIterator;

    public ArtificialMultiSampleReadStream( Collection<ArtificialSingleSampleReadStream> perSampleArtificialReadStreams ) {
        if ( perSampleArtificialReadStreams == null || perSampleArtificialReadStreams.isEmpty() ) {
            throw new ReviewedStingException("Can't create an ArtificialMultiSampleReadStream out of 0 ArtificialSingleSampleReadStreams");
        }

        this.perSampleArtificialReadStreams = perSampleArtificialReadStreams;
    }

    public Iterator<SAMRecord> iterator() {
        // lazy initialization to prevent reads from being created until they're needed
        initialize();

        return mergingIterator;
    }

    public StingSAMIterator getStingSAMIterator() {
        // lazy initialization to prevent reads from being created until they're needed
        initialize();

        return StingSAMIteratorAdapter.adapt(mergingIterator);
    }

    private void initialize() {
        Collection<SAMFileReader> perSampleSAMReaders = new ArrayList<SAMFileReader>(perSampleArtificialReadStreams.size());
        Collection<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(perSampleArtificialReadStreams.size());

        for ( ArtificialSingleSampleReadStream readStream : perSampleArtificialReadStreams ) {
            Collection<SAMRecord> thisStreamReads = readStream.makeReads();

            SAMFileReader reader = new ArtificialSAMFileReader(readStream.getHeader(),
                                                               thisStreamReads.toArray(new SAMRecord[thisStreamReads.size()]));
            perSampleSAMReaders.add(reader);
            headers.add(reader.getFileHeader());
        }

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, true);
        mergingIterator = new MergingSamRecordIterator(headerMerger, perSampleSAMReaders, true);
    }
}

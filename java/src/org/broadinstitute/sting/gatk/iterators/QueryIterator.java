package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordComparator;
import net.sf.samtools.SAMFileReader;

import java.util.Iterator;

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
 * @author aaron
 *
 * This interface indicates that the iterator is query able
 */
public interface QueryIterator  extends StingSAMIterator, PeekingIterator<SAMRecord> {

    /**
     * The required methods to query able
     **/
    public void queryOverlapping(final String contig, final int start, final int stop);
    public void query(final String contig, final int start, final int stop, final boolean contained);
    public void queryUnmappedReads();
    public void queryContained(final String contig, final int start, final int stop);
}

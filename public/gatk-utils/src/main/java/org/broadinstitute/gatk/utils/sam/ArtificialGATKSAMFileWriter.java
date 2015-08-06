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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;

import java.util.ArrayList;
import java.util.List;


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
 *         <p/>
 *         Class ArtificialGATKSAMFileWriter
 *         <p/>
 * generates a fake samwriter, that you can get the output reads
 * from when you're done.  
 */
public class ArtificialGATKSAMFileWriter implements GATKSAMFileWriter {

    // are we closed
    private boolean closed = false;

    // the SAMRecords we've added to this writer
    List<SAMRecord> records = new ArrayList<SAMRecord>();

    public void addAlignment( SAMRecord alignment ) {
        records.add(alignment);
    }

    public SAMFileHeader getFileHeader() {
        if (records.size() > 0) {
            return records.get(0).getHeader();
        }
        return null;
    }

    /** not much to do when we're fake */
    public void close() {
        closed = true;
    }

    /**
     * are we closed?
     *
     * @return true if we're closed
     */
    public boolean isClosed() {
        return closed;
    }

    /**
     * get the records we've seen
     * @return
     */
    public List<SAMRecord> getRecords() {
        return records;
    }

    @Override
    public void writeHeader(SAMFileHeader header) {
    }

    @Override
    public void setPresorted(boolean presorted) {
    }

    @Override
    public void setMaxRecordsInRam(int maxRecordsInRam) {
    }

    /**
     * @throws java.lang.UnsupportedOperationException No progress logging in this implementation.
     */
    @Override
    public void setProgressLogger(final ProgressLoggerInterface logger) {
        throw new UnsupportedOperationException("Progress logging not supported");
    }
}

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
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;

/**
 * XXX
 */
public class SimplifyingSAMFileWriter implements SAMFileWriter {
    final SAMFileWriter dest;

    public SimplifyingSAMFileWriter(final SAMFileWriter finalDestination) {
        this.dest = finalDestination;
    }

    public void addAlignment( SAMRecord read ) {
        if ( keepRead(read) ) {
            dest.addAlignment(simplifyRead(read));

        }
    }

    /**
     * Retrieves the header to use when creating the new SAM file.
     * @return header to use when creating the new SAM file.
     */
    public SAMFileHeader getFileHeader() {
        return dest.getFileHeader();
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        dest.close();
    }


    public static final boolean keepRead(SAMRecord read) {
        return ! excludeRead(read);
    }

    public static final boolean excludeRead(SAMRecord read) {
        return read.getReadUnmappedFlag() || read.getReadFailsVendorQualityCheckFlag() || read.getDuplicateReadFlag() || read.getNotPrimaryAlignmentFlag();
    }

    public static final SAMRecord simplifyRead(SAMRecord read) {
        // the only attribute we keep is the RG
        Object rg = read.getAttribute("RG");
        read.clearAttributes();
        read.setAttribute("RG", rg);
        return read;
    }

    @Override
    public void setProgressLogger(final ProgressLoggerInterface logger) {
        dest.setProgressLogger(logger);
    }
}

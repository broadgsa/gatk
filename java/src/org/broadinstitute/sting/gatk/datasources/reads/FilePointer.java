/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.SAMFileSpan;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
class FilePointer {
    protected final Map<SAMReaderID,SAMFileSpan> fileSpans = new HashMap<SAMReaderID,SAMFileSpan>();
    protected final String referenceSequence;
    protected final BAMOverlap overlap;
    protected final List<GenomeLoc> locations;

    /**
     * Does this file pointer point into an unmapped region?
     */
    protected final boolean isRegionUnmapped;

    public FilePointer(final GenomeLoc location) {
        this.referenceSequence = location.getContig();
        this.overlap = null;
        this.locations = Collections.singletonList(location);
        this.isRegionUnmapped = GenomeLoc.isUnmapped(location);
    }

    public FilePointer(final String referenceSequence,final BAMOverlap overlap) {
        this.referenceSequence = referenceSequence;
        this.overlap = overlap;
        this.locations = new ArrayList<GenomeLoc>();
        this.isRegionUnmapped = false;
    }

    public void addLocation(GenomeLoc location) {
        locations.add(location);
    }

    public void addFileSpans(SAMReaderID id, SAMFileSpan fileSpan) {
        this.fileSpans.put(id,fileSpan);
    }
}

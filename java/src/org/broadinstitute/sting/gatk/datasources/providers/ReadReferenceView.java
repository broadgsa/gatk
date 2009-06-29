package org.broadinstitute.sting.gatk.datasources.providers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequence;
import org.broadinstitute.sting.utils.Utils;
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
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:36:14 PM
 *
 */

/** Provides access to the reference over a single read. */

public class ReadReferenceView extends ReferenceView {
    /**
     * Create a view of the reference with respect to a single read.
     *
     * @param provider
     */
    public ReadReferenceView( ShardDataProvider provider ) {
        super(provider);
    }

    /**
     * Gets the bases of the reference that are aligned to the given read.
     *
     * @param read the read for which to extract reference information.
     *
     * @return The bases corresponding to this read, or null if the read is unmapped.
     *         If the alignment goes off the end of the contig, return just the portion
     *         mapped to the reference, followed by X's coresponding to the rest of the read.
     *         This indicates that the rest lies off the end of the contig.
     */
    public char[] getReferenceBases( SAMRecord read ) {
        if (read.getReadUnmappedFlag())
            return null;

        String contig = read.getReferenceName();
        int start = read.getAlignmentStart();
        int stop = read.getAlignmentEnd();

        SAMSequenceRecord sequenceRecord = reference.getSequenceDictionary().getSequence(contig);
        if (stop > sequenceRecord.getSequenceLength())
            stop = sequenceRecord.getSequenceLength();

        ReferenceSequence alignmentToReference = reference.getSubsequenceAt(contig, start, stop);
        return ( StringUtil.bytesToString(alignmentToReference.getBases()) + Utils.dupString('X', read.getAlignmentEnd() - stop) ).toCharArray();
    }

}

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

package org.broadinstitute.sting.gatk.contexts;

import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import net.sf.samtools.SAMRecord;
import java.util.*;


public class VariantContext {
    private RefMetaDataTracker tracker;
    private ReferenceContext ref;
    private AlignmentContext context;
    private rodVariants variant;
    private AlignmentContext Q0freeContext = null;

    public VariantContext(RefMetaDataTracker tracker, ReferenceContext ref,
                          AlignmentContext context, rodVariants variant) {
        this.tracker = tracker;
        this.ref = ref;
        this.context = context;
        this.variant = variant;
    }

    public RefMetaDataTracker getTracker() { return tracker; }
    public ReferenceContext getReferenceContext() { return ref; }
    public rodVariants getVariant() { return variant; }
    public AlignmentContext getAlignmentContext() { return getAlignmentContext(false); }
    public AlignmentContext getAlignmentContext(boolean useMQ0Reads) {
        return (useMQ0Reads ? context : getQ0freeContext());
    }

    private AlignmentContext getQ0freeContext() {
        if ( Q0freeContext == null ) {
            // set up the variables
            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();
            Iterator<SAMRecord> readIter = reads.iterator();
            Iterator<Integer> offsetIter = offsets.iterator();
            ArrayList<SAMRecord> Q0freeReads = new ArrayList<SAMRecord>();
            ArrayList<Integer> Q0freeOffsets = new ArrayList<Integer>();

            // copy over good reads/offsets
            while ( readIter.hasNext() ) {
                SAMRecord read = readIter.next();
                Integer offset = offsetIter.next();
                if ( read.getMappingQuality() > 0 ) {
                    Q0freeReads.add(read);
                    Q0freeOffsets.add(offset);
                }
            }

            Q0freeContext = new AlignmentContext(context.getLocation(), Q0freeReads, Q0freeOffsets);
        }

        return Q0freeContext;
    }
}
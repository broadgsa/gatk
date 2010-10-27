/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.*;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

/**
 * Left aligns indels in reads.
 */
public class LeftAlignIndels extends ReadWalker<Integer, Integer> {

    @Output(required=false, doc="Output bam")
    protected StingSAMFileWriter writer = null;

    @Argument(fullName="maxReadsInRam", shortName="maxInRam", doc="max reads allowed to be kept in memory at a time by the SAMFileWriter. "+
                "If too low, the tool may run out of system file descriptors needed to perform sorting; if too high, the tool may run out of memory.", required=false)
    protected int MAX_RECORDS_IN_RAM = 500000;

    public void initialize() {
        // set up the output writer
        if ( writer != null )
            writer.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);
    }

    private void emit(final SAMRecord read) {
        if ( writer != null )
            writer.addAlignment(read);
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        // we can not deal with screwy records
        if ( read.getCigar().numCigarElements() == 0 ) {
            emit(read);
            return 0;
        }

        // move existing indels (for 1 indel reads only) to leftmost position within identical sequence
        int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
        if ( numBlocks == 2 ) {
            Cigar newCigar = AlignmentUtils.leftAlignIndel(IndelRealigner.unclipCigar(read.getCigar()), ref.getBases(), read.getReadBases(), 0, 0);
            newCigar = IndelRealigner.reclipCigar(newCigar, read);
            read.setCigar(newCigar);
        }

        emit(read);
        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {}
}
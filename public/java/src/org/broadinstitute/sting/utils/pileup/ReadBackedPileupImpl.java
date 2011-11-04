/*
 * Copyright (c) 2010, The Broad Institute
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
package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.List;
import java.util.Map;

public class ReadBackedPileupImpl extends AbstractReadBackedPileup<ReadBackedPileupImpl,PileupElement> implements ReadBackedPileup {

    public ReadBackedPileupImpl(GenomeLoc loc) {
        super(loc);
    }

    public ReadBackedPileupImpl(GenomeLoc loc, List<GATKSAMRecord> reads, List<Integer> offsets ) {
        super(loc,reads,offsets);
    }

    public ReadBackedPileupImpl(GenomeLoc loc, List<GATKSAMRecord> reads, int offset ) {
        super(loc,reads,offset);
    }

    public ReadBackedPileupImpl(GenomeLoc loc, List<PileupElement> pileupElements) {
        super(loc,pileupElements);
    }

    public ReadBackedPileupImpl(GenomeLoc loc, Map<String,ReadBackedPileupImpl> pileupElementsBySample) {
        super(loc,pileupElementsBySample);
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc
     * @param pileup
     */
    public ReadBackedPileupImpl(GenomeLoc loc, List<PileupElement> pileup, int size, int nDeletions, int nMQ0Reads) {
        super(loc,pileup,size,nDeletions,nMQ0Reads);
    }

    protected ReadBackedPileupImpl(GenomeLoc loc, PileupElementTracker<PileupElement> tracker) {
        super(loc,tracker);
    }

    @Override
    protected ReadBackedPileupImpl createNewPileup(GenomeLoc loc, PileupElementTracker<PileupElement> tracker) {
        return new ReadBackedPileupImpl(loc,tracker);
    }

    @Override
    protected PileupElement createNewPileupElement(GATKSAMRecord read, int offset) {
        return new PileupElement(read,offset);
    }
}

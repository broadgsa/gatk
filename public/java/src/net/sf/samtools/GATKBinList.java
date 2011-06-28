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

package net.sf.samtools;

import java.util.BitSet;

/**
 * A temporary solution to work around Java access rights issues:
 * override chunk and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKBinList extends BinList {
    /**
     * Create a new BinList over sequenceCount sequences, consisting of the given bins.
     * @param referenceSequence Reference sequence to which these bins are relevant.
     * @param bins The given bins to include.
     */
    public GATKBinList(final int referenceSequence, final BitSet bins) {
        super(referenceSequence,bins);
    }

    /**
     * Retrieves the bins stored in this list.
     * @return A bitset where a bin is present in the list if the bit is true.
     */
    public BitSet getBins() {
        return super.getBins();
    }
}

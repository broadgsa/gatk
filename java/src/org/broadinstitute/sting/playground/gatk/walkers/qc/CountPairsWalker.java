/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.qc;

import org.broadinstitute.sting.gatk.walkers.ReadPairWalker;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import net.sf.samtools.SAMRecord;

import java.util.Collection;
import java.util.List;

/**
 * Counts the number of read pairs encountered in a file sorted in
 * query name order.  Breaks counts down by total pairs and number
 * of paired reads.
 *
 * @author mhanna
 * @version 0.1
 */
public class CountPairsWalker extends ReadPairWalker<Integer,Long> {
    /**
     * How many reads are the first in a pair, based on flag 0x0040 from the SAM spec.
     */
    private long firstOfPair = 0;

    /**
     * How many reads are the second in a pair, based on flag 0x0080 from the SAM spec.
     */
    private long secondOfPair = 0;

    /**
     * A breakdown of the total number of reads seen with exactly the same read name.
     */
    private List<Long> pairCountsByType = new ExpandingArrayList<Long>();

    /**
     * Maps a read pair to a given reduce of type MapType.  Semantics determined by subclasser.
     * @param reads Collection of reads having the same name.
     * @return Semantics defined by implementer.
     */
    @Override
    public Integer map(Collection<SAMRecord> reads) {
        if(pairCountsByType.get(reads.size()) != null)
            pairCountsByType.set(reads.size(),pairCountsByType.get(reads.size())+1);
        else
            pairCountsByType.set(reads.size(),1L);

        for(SAMRecord read: reads) {
            if(read.getFirstOfPairFlag()) firstOfPair++;
            if(read.getSecondOfPairFlag()) secondOfPair++;
        }

        return 1;
    }

    /**
     * No pairs at the beginning of a traversal.
     * @return 0 always.
     */
    @Override
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Combine number of pairs seen in this iteration (always 1) with total number of pairs
     * seen in previous iterations. 
     * @param value Pairs in this iteration (1), from the map function.
     * @param sum Count of all pairs in prior iterations.
     * @return All pairs encountered in previous iterations + all pairs encountered in this iteration (sum + 1).
     */
    @Override
    public Long reduce(Integer value, Long sum) {
        return value + sum;    
    }

    /**
     * Print summary statistics over the entire traversal.
     * @param sum A count of all read pairs viewed.
     */
    @Override
    public void onTraversalDone(Long sum) {
        out.printf("Total number of pairs               : %d%n",sum);
        out.printf("Total number of first reads in pair : %d%n",firstOfPair);
        out.printf("Total number of second reads in pair: %d%n",secondOfPair);
        for(int i = 1; i < pairCountsByType.size(); i++) {
            if(pairCountsByType.get(i) == null)
                continue;
            out.printf("Pairs of size %d: %d%n",i,pairCountsByType.get(i));
        }
    }

}

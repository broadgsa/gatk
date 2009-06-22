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

package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * a class to store the duplicates information we pass around
 */
class DuplicateCount  {
        public int count = 0;
        public int undupDepth = 0;
        public int depth = 0;
    }

/**
 * 
 * @author aaron 
 * 
 * Class CountDuplicatesWalker
 *
 * Count the number of duplicates, and the average depth of duplicates at all positions.
 */
public class CountDuplicatesWalker extends DuplicateWalker<DuplicateCount, DuplicateCount> {

    /**
     * the map function, conforming to the duplicates interface
     * @param loc the genomic location
     * @param refBases the reference bases that cover this position, which we turn off
     * @param context the LocusContext, containing all the reads overlapping this region
     * @param uniqueReads all the unique reads, bundled together
     * @param duplicateReads all the duplicate reads
     * @return a DuplicateCount object, with the appropriate stats
     */
    public DuplicateCount map(GenomeLoc loc, byte[] refBases, LocusContext context, List<SAMRecord> uniqueReads, List<SAMRecord> duplicateReads) {
        DuplicateCount dup = new DuplicateCount();
        dup.count = 1;
        dup.undupDepth = uniqueReads.size();
        dup.depth = duplicateReads.size();
        return dup;
    }

    /**
     * setup our walker.  In this case, new a DuplicateCount object and return it
     * @return the object holding the counts of the duplicates
     */
    public DuplicateCount reduceInit() {
        return new DuplicateCount();
    }

    /**
     * the reduce step.  This function combines the DuplicateCount objects, and updates the running depth average
     * @param value the new DuplicateCount
     * @param sum the running sum DuplicateCount
     * @return a new DuplicateCount with the updated sums
     */
    public DuplicateCount reduce(DuplicateCount value, DuplicateCount sum) {
        DuplicateCount dup = new DuplicateCount();
        dup.count = sum.count + value.count;
        dup.depth = value.depth + sum.depth;
        dup.undupDepth = value.undupDepth + sum.undupDepth;
        return dup;
    }

    /**
     * when we're done, print out the collected stats
     * @param result the result of the traversal engine, to be printed out
     */
    public void onTraversalDone(DuplicateCount result) {
        out.println("[REDUCE RESULT] Traversal result is: ");
        out.println("duplicate count = " + result.count);
        out.println("average depth = " + (double)result.depth / (double)result.count);
        out.println("duplicates seen = " + result.depth);
        out.println("average unique depth = " + (double)result.undupDepth / (double)result.undupDepth);
        out.println("unique read count = " + result.undupDepth);

    }
}

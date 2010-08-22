/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.io.PrintStream;

/**
 * a class to store the traversal information we pass around
 */
class DuplicateCount  {
    public int count = 0;       // the count of sites we were given
    public int nUniqueMolecules = 0;  // the unique read count
    public int nDuplicatedMolecules = 0;  // the unique read count
    public int depth = 0;       // the dupplicate read depth
}

/**
 * Count the number of unique reads, duplicates, and the average depth of unique reads and duplicates at all positions.
 * @author mark DePristo
 */
public class CountDuplicatesWalker extends DuplicateWalker<DuplicateCount, DuplicateCount> {
    @Output
    PrintStream out;

    @Argument(fullName="quietLocus", required=false, doc="If true, per locus information isn't printed")
    public boolean quiet = false;

    /**
     * the map function, conforming to the duplicates interface
     * @param loc the genomic location
     * @param context the AlignmentContext, containing all the reads overlapping this region
     * @param readSets all the duplicate reads
     * @return a DuplicateCount object, with the appropriate stats
     */
    public DuplicateCount map(GenomeLoc loc, AlignmentContext context, Set<List<SAMRecord>> readSets ) {
        if ( ! quiet ) out.printf("%s with %d read sets => ", loc, readSets.size());

        DuplicateCount dup = new DuplicateCount();
        dup.depth = 0;
        for ( List<SAMRecord> reads : readSets) {
            List<String> names = new ArrayList<String>();
            for ( SAMRecord read : reads ) {
                names.add(read.getReadName());
            }
            if ( ! quiet ) out.printf("%d reads [%s] ", reads.size(), Utils.join(",", names));
            dup.depth += reads.size();
            dup.nDuplicatedMolecules += reads.size() > 1 ? 1 : 0;   // if there's more than 1 read per set, we're a duplicated reads
        }
        if ( ! quiet ) out.printf("%n");

        dup.count = 1;
        dup.nUniqueMolecules = readSets.size();
        return dup;
    }

    public boolean mapAtLociWithoutDuplicates() { return true; }

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
        dup.nDuplicatedMolecules = value.nDuplicatedMolecules + sum.nDuplicatedMolecules;
        dup.nUniqueMolecules = value.nUniqueMolecules + sum.nUniqueMolecules;
        return dup;
    }

    /**
     * when we're done, print out the collected stats
     * @param result the result of the traversal engine, to be printed out
     */
    public void onTraversalDone(DuplicateCount result) {
        out.println("[REDUCE RESULT] Traversal result is: ");
        out.println("traversal iterations = " + result.count);
        out.printf("average depth = %.2f%n", (double)result.depth / (double)result.count);
        out.println("unique molecules seen = " + result.nUniqueMolecules);
        out.println("duplicated molecules seen = " + result.nDuplicatedMolecules);
        out.printf("percent duplicated = %.2f%%%n", result.nDuplicatedMolecules / (double)result.nUniqueMolecules * 100);
        out.printf("average unique read depth = %.2f%n", (double)result.nUniqueMolecules / (double)result.count);
    }
}

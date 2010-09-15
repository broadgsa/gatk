package org.broadinstitute.sting.playground.sample;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Walks over the input data set, calculating the total number of covered loci for diagnostic purposes.
 * Simplest example of a locus walker.
 */
public class CountLociByPopulationWalker extends LocusWalker<Integer, Long> implements TreeReducible<Long> {
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        HashMap<String, Integer> count = new HashMap<String, Integer>();

        ArrayList<SAMRecord> reads = (ArrayList) context.getBasePileup().getReads();

        for (SAMRecord read : reads) {
            String population = getToolkit().getSampleByRead(read).getPopulation();
            if (!count.containsKey(population)) {
                count.put(population, 1);
            }
            count.put(population, count.get(population) + 1);
        }

        System.out.println("\n\n\n*****  LOCUS: " + ref.toString() + "  *****");
        for (String population : count.keySet()) {
            System.out.println(String.format("%s | %d\n", population, count.get(population)));
        }

        return 1;
    }

    public Long reduceInit() { return 0l; }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }

    /**
     * Reduces two subtrees together.  In this case, the implementation of the tree reduce
     * is exactly the same as the implementation of the single reduce.
     */
    public Long treeReduce(Long lhs, Long rhs) {
        return lhs + rhs;
    }
}
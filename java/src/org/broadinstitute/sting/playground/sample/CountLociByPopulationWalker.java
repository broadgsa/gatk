package org.broadinstitute.sting.playground.sample;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Extends locus walker to print how many reads there are at each locus, by population
 */
public class CountLociByPopulationWalker extends LocusWalker<Integer, Long> implements TreeReducible<Long> {
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        // in this HashMap, we'll keep count of how many
        HashMap<String, Integer> count = new HashMap<String, Integer>();

        ArrayList<SAMRecord> reads = (ArrayList) context.getBasePileup().getReads();

        for (SAMRecord read : reads) {

            // get the sample
            Sample sample = getToolkit().getSampleByRead(read);
            if (sample == null)
                return 1;

            if (!count.containsKey(sample.getPopulation())) {
                count.put(sample.getPopulation(), 1);
            }
            count.put(sample.getPopulation(), count.get(sample.getPopulation()) + 1);
        }

        System.out.println("\n\n\n*****  LOCUS: " + ref.getLocus().toString() + "  *****");
        for (String population : count.keySet()) {
            System.out.println(String.format("%s | %d", population, count.get(population)));
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
package org.broadinstitute.sting.playground.contexts;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.Pair;

import java.util.List;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 9, 2009
 * Time: 10:43:23 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FilteredAlignmentContext extends AlignmentContext{

    public FilteredAlignmentContext() { /* super method is called */ }

    /* A partitioned alignment context must have a constructor
     * method that generates the object from another alignment
     * context
     */

    public FilteredAlignmentContext(AlignmentContext context) {
        Pair<List<SAMRecord>, List<Integer>> partitionedReads = filter(context);
        this.reads = partitionedReads.getFirst();
        this.offsets = partitionedReads.getSecond();
        this.loc = context.getLocation();

    }

    /*
     * A new partitioned alignment object need only specify how the reads from an Alignmentcontext
     * are to be partitioned, and return the new partition in a pair.
     * @Param: context - an alignment context containing reads to be partitioned
     */
    public abstract Pair<List<SAMRecord>, List<Integer>> filter(AlignmentContext context);

}

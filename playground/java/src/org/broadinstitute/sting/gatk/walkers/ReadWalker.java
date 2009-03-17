package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ReadWalker<MapType, ReduceType> extends Walker {
    boolean requiresOrderedReads();
    
    // Do we actually want to operate on the context?
    boolean filter(LocusContext context, SAMRecord read);

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    MapType map(LocusContext context, SAMRecord read);

    // Given result of map function
    ReduceType reduceInit();
    ReduceType reduce(MapType value, ReduceType sum);
}

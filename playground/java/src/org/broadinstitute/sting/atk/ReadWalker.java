package org.broadinstitute.sting.atk;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.atk.LocusContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ReadWalker<MapType, ReduceType> {
    void initialize();
    public String walkerType();
    
    // Do we actually want to operate on the context?
    boolean filter(LocusContext context, SAMRecord read);

    // Map over the org.broadinstitute.sting.atk.LocusContext
    MapType map(LocusContext context, SAMRecord read);

    // Given result of map function
    ReduceType reduceInit();
    ReduceType reduce(MapType value, ReduceType sum);

    void onTraversalDone();
}

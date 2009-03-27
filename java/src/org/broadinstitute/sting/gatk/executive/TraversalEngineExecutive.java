package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 27, 2009
 * Time: 10:00:05 AM
 * To change this template use File | Settings | File Templates.
 */
public interface TraversalEngineExecutive<ReduceType> {
    
    public ReduceType processIntervals(List<GenomeLoc> locations);
}

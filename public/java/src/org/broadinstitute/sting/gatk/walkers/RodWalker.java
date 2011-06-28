package org.broadinstitute.sting.gatk.walkers;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE_ORDERED_DATA)
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_ORDERED_DATA})
@Allows({DataSource.REFERENCE, DataSource.REFERENCE_ORDERED_DATA})
public abstract class RodWalker<MapType, ReduceType> extends LocusWalker<MapType, ReduceType> {
}
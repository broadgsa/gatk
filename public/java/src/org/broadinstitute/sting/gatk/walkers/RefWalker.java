package org.broadinstitute.sting.gatk.walkers;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@Allows(DataSource.REFERENCE)
public abstract class RefWalker<MapType, ReduceType> extends LocusWalker<MapType, ReduceType> {
}
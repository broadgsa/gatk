package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.LocusContext;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RefWalker<MapType, ReduceType> extends LocusWalker<MapType, ReduceType> {
    public boolean requiresReads()     { return false; }
    public boolean cannotHandleReads() { return true; }
}
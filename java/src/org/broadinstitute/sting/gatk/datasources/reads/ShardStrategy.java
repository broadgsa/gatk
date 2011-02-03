package org.broadinstitute.sting.gatk.datasources.reads;

import java.util.Iterator;
/**
 *
 * User: aaron
 * Date: Apr 10, 2009
 * Time: 4:55:37 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */

/**
 * @author aaron
 * @version 1.0
 * @date Apr 10, 2009
 * <p/>
 * Interface ShardStrategy
 * <p/>
 * The base interface for the sharding strategy; before we had a base abstract
 * class, but not this will be an interface to accomidate read based sharding
 */
public interface ShardStrategy extends Iterator<Shard>, Iterable<Shard> {
}

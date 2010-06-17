package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.Serializable;
import java.util.List;
/**
 *
 * User: aaron
 * Date: Apr 10, 2009
 * Time: 5:00:27 PM
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
 * Interface Shard
 * <p/>
 * The base interface for shards.
 */
public interface Shard extends Serializable {
    enum ShardType {
        READ, LOCUS
    }

    /** @return the genome location represented by this shard */
    public List<GenomeLoc> getGenomeLocs();

    /**
     * what kind of shard do we return
     * @return ShardType, indicating the type
     */
    public ShardType getShardType();
}

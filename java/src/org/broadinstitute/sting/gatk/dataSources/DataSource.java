package org.broadinstitute.sting.gatk.dataSources;

import org.broadinstitute.sting.gatk.dataSources.shards.DataShard;

/**
 * User: aaron
 * Date: Mar 25, 2009
 * Time: 6:20:00 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public interface DataSource {

    public DataShard toChunk(int chunkCount);
}

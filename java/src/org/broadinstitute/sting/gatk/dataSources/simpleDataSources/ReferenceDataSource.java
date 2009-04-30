package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.BoundedReferenceIterator;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 3:55:21 PM
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
 * @date Apr 6, 2009
 * <p/>
 * Class ReferenceDataSource
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class ReferenceDataSource implements SimpleDataSource {

    final protected IndexedFastaSequenceFile refFile;

    /**
     * Query the data source for a region of interest, specified by the genome location.
     * The iterator will generate successive calls
     *
     * @param shard the genome location to extract data for
     * @return an iterator of the appropriate type, that is limited by the region
     */
    public BoundedReferenceIterator seek(Shard shard) {
        if (shard.getShardType() == Shard.ShardType.LOCUS) {
            BoundedReferenceIterator ret = new BoundedReferenceIterator(refFile, shard.getGenomeLoc());
            return ret;
        } else {
            throw new StingException("ReferenceDataSource can only take LocusShards");
        }

    }

    public ReferenceDataSource(String refFileName) throws SimpleDataSourceLoadException {
        if (refFileName == null) {
            throw new SimpleDataSourceLoadException("ReferenceDataSource: refFileName passed in is null");
        }
        File infile = new File(refFileName);
        if (!infile.canRead()) {
            throw new SimpleDataSourceLoadException("ReferenceDataSource: Unable to load file: " + refFileName);
        }
        try {
            refFile = new IndexedFastaSequenceFile(new File(refFileName));
        }
        catch( FileNotFoundException ex ) {
            throw new SimpleDataSourceLoadException( "Unable to find reference file", ex );
        }
    }
}

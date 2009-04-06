package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.FastaSequenceFile2;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;

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

    final protected FastaSequenceFile2 refFile;
    final protected ReferenceIterator refIter;

    /**
     * Query the data source for a region of interest, specified by the genome location.
     * The iterator will generate successive calls
     *
     * @param location the genome location to extract data for
     * @return an iterator of the appropriate type, that is limited by the region
     */
    public ReferenceIterator seek(GenomeLoc location) {
        ReferenceIterator refSite = refIter.seekForward(location);
        return refSite;
    }

    /**
     * Constructor - ReferenceDataSource
     *
     * @param refFileName the reference file
     * @throws SimpleDataSourceLoadException
     */
    public ReferenceDataSource(String refFileName) throws SimpleDataSourceLoadException {
        if (refFileName == null) {
            throw new SimpleDataSourceLoadException("ReferenceDataSource: refFileName passed in is null");
        }
        File infile = new File(refFileName);
        if (!infile.canRead()) {
            throw new SimpleDataSourceLoadException("ReferenceDataSource: Unable to load file: " + refFileName);
        }
        refFile = new FastaSequenceFile2(new File(refFileName));
        refIter = new ReferenceIterator(this.refFile);

    }
}

package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import org.broadinstitute.sting.gatk.refdata.HapMapAlleleFrequenciesROD;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.*;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 4:33:10 PM
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
 * Class ReferenceMetaDataSource
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class ReferenceMetaDataSource implements SimpleDataSource {

    // our enumerated types
    public enum RODTYPE {
        DBSNP, HAPMAP
    }

    // these could go on the stack, but a heap copy isn't too bad
    private List<ReferenceOrderedDatum> myData = null;
    private List<ReferenceOrderedData<? extends ReferenceOrderedDatum>.RODIterator> rodIters = null;
    private List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods = null;

    /**
     * Prepare the list of reference ordered data iterators for each of the rods
     *
     * @return A list of ROD iterators for getting data from each ROD
     */
    protected List<ReferenceOrderedData<? extends ReferenceOrderedDatum>.RODIterator> initializeRODs() {
        // set up reference ordered data
        rodIters = new ArrayList<ReferenceOrderedData<? extends ReferenceOrderedDatum>.RODIterator>();
        for (ReferenceOrderedData<? extends ReferenceOrderedDatum> data : rods) {
            rodIters.add(data.iterator());
        }
        return rodIters;
    }

    /**
     * Builds a list of the reference ordered datum at loc from each of the iterators.  This function
     * assumes you are accessing the data in order.  You can't use this function for random access.  Each
     * successive call moves you along the file, consuming all data before loc.
     *
     * @param rodIters Iterators to access the RODs
     * @param loc      The location to get the rods at
     * @return A list of ReferenceOrderDatum at loc.  ROD without a datum at loc will be null in the list
     */
    protected List<ReferenceOrderedDatum> getReferenceOrderedDataAtLocus(List<ReferenceOrderedData<? extends ReferenceOrderedDatum>.RODIterator> rodIters,
                                                                         final GenomeLoc loc) {
        List<ReferenceOrderedDatum> data = new ArrayList<ReferenceOrderedDatum>();
        for (ReferenceOrderedData<? extends ReferenceOrderedDatum>.RODIterator iter : rodIters) {
            data.add(iter.seekForward(loc));
        }
        return data;
    }

    /**
     * Query the data source for a region of interest, specified by the genome location.
     * The iterator will generate successive calls
     *
     * @param location the genome location to extract data for
     * @return an iterator of the appropriate type, that is limited by the region
     */
    public Iterator<ReferenceOrderedDatum> seek(GenomeLoc location) {
        myData = getReferenceOrderedDataAtLocus(rodIters, location);
        return myData.iterator();
    }

    public ReferenceMetaDataSource(HashMap<String, RODTYPE> files) {

        // setup a rod list
        List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods = new ArrayList<ReferenceOrderedData<? extends ReferenceOrderedDatum>>();

        // cycle through the passed in rod's

        Set<String> fileNames = files.keySet();
        for (String file : fileNames) {
            switch (files.get(file)) {

                case DBSNP: {
                    ReferenceOrderedData<rodDbSNP> dbsnp = new ReferenceOrderedData<rodDbSNP>(new File(file), rodDbSNP.class);
                    //dbsnp.testMe();
                    rods.add(dbsnp); // { gff, dbsnp };
                }
                case HAPMAP: {
                    ReferenceOrderedData<HapMapAlleleFrequenciesROD> hapmap = new ReferenceOrderedData<HapMapAlleleFrequenciesROD>(new File(file), HapMapAlleleFrequenciesROD.class);
                    //dbsnp.testMe();
                    rods.add(hapmap); // { gff, dbsnp };
                }
            }
        }
    }
}

package org.broadinstitute.sting.gatk.refdata.utils;

import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;


/**
 * 
 * @author aaron 
 * 
 * Class RMDIntervalGenerator
 *
 * Creates an interval list, given an RMDTrack
 */
public class RMDIntervalGenerator {
    public ReferenceOrderedDataSource dataSource;

    /**
     * create a interval representation of a ROD track
     * @param dataSource the track
     */
    public RMDIntervalGenerator(ReferenceOrderedDataSource dataSource) {
        if (dataSource == null) throw new IllegalArgumentException("Data source cannot be null");
        this.dataSource = dataSource;
    }

    /**
     * create a genome location list from the interval track
     * @return a list of genome locations
     */
    public List<GenomeLoc> toGenomeLocList() {
        Iterator<RODRecordList> iter = dataSource.seek((GenomeLoc)null);
        List<GenomeLoc> locations = new ArrayList<GenomeLoc>();
        while (iter.hasNext()) {
            RODRecordList feature = iter.next();
            GenomeLoc loc = feature.getLocation();
            if (loc != null) locations.add(loc);            
        }
        return locations;
    }

    /**
     * return a map of reference meta data track names to RODS
     * @param sources the reference ordered data sources to get the names from
     * @return a map of reference meta data names to RODS
     */
    public static Map<String,ReferenceOrderedDataSource> getRMDTrackNames(List<ReferenceOrderedDataSource> sources) {
        // get a list of the current rod names we're working with
        Map<String,ReferenceOrderedDataSource> rodNames = new HashMap<String,ReferenceOrderedDataSource>();
        for (ReferenceOrderedDataSource rod : sources) {
            rodNames.put(rod.getName(),rod);
        }
        return rodNames;
    }
}

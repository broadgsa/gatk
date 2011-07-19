package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 21, 2009
 * Time: 2:54:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class PlatformUnitFilter extends ReadFilter {
    // a hack: use static in order to be able to fill it with the data from command line at runtime
    static private Set<String> blackListedLanes = new HashSet<String>();

    public boolean filterOut(SAMRecord samRecord) {

        if ( blackListedLanes.size() == 0 ) return false; // no filters set, nothing to do

        Object pu_attr = samRecord.getAttribute("PU");

        if ( pu_attr == null ) {
            // no platform unit in the record, go get from read group
            SAMReadGroupRecord rgr = samRecord.getReadGroup();
            if ( rgr == null ) throw new UserException.MalformedBAM(samRecord, "Read " + samRecord.getReadName() +" has NO associated read group record");
            pu_attr = rgr.getAttribute("PU") ;
        }
        if ( pu_attr == null ) return false; // could not get PU, forget about the filtering...
        return blackListedLanes.contains((String)pu_attr);
    }

    /**
     * The argument is interpreted as a comma-separated list of lanes (platform units) to be filtered
     * out. All the specified names will be registered with the filter and filterOut(r) for any SAMRecord r
     * belonging to one of the specified lanes will thereafter return true.
     * The names can be surrounded by additional spaces, the latters will be trimmed by this method.
     * This method can be called multiple times to add more lanes. Re-registering the same lane again is safe.
     * @param arg
     */
    public static void setBlackListedLanes(String arg) {
        String[] lanes = arg.split(",");
        for ( int i = 0; i < lanes.length ; i++ ) {
            blackListedLanes.add(lanes[i].trim());
        }
    }

    /**
     * Adds a single name of a lane (platform unit) to be filtered out by this filter. The name can be surrounded
     * by spaces, the latters will be trimmed out. This method can be called multiple times to add more lanes.
     * Re-registering the same lane again is safe.
     * @param arg
     */
    public static void addBlackListedLane(String arg) {
        blackListedLanes.add(arg.trim());
    }

}

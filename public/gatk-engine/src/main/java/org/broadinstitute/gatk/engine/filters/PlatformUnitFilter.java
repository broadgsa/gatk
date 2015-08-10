/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.filters;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.HashSet;
import java.util.Set;

/**
 * Filter out reads with blacklisted platform unit tags
 *
 * <p>This filter is useful for running on only a subset of the data as identified by a read group property.
 * In the case of the platform unit filter, the goal is usually to blacklist certain runs if we know there was a problem with
 * a particular sequencing machine.</p>
 *
 * @author asivache
 * @since Sep 21, 2009
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

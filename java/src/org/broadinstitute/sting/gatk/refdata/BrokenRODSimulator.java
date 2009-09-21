package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;

import java.util.Map;
import java.util.HashMap;

/**
 * This is a temporary solution that keeps integration tests passing until they are fixed permanently.
 * The new ROD iterator system fixes a few issues present in the previous version, and as the result applications
 * see somewhat different sets of RODs (notably, in ubiquitous rodDbSNP). This class takes the results returned
 * by the new ROD system and simulates the results that would be returned by the old one. Everytime this class is used,
 * it's an indication of the urgent need to get rid of it and fix the integration test!!!
 *
 *
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 21, 2009
 * Time: 11:15:34 AM
 * To change this template use File | Settings | File Templates.
 */
public class BrokenRODSimulator {
    // multiple instances can access the sam tracker (and thus underlying iterator) from different
    // places in the code. By making the mapping static we simulate this paradigm of "accessing the same RODIterator"
    // through the tracke objects passed around at will.
    private static Map<String,GenomeLoc> last_intervals = new HashMap<String,GenomeLoc>();
    private static Map<String,ReferenceOrderedDatum> last_rods = new HashMap<String,ReferenceOrderedDatum>();

    public BrokenRODSimulator() {
//        last_interval = GenomeLocParser.createGenomeLoc(0,1,1);
    }

    public static void attach(String name) {
        if ( last_intervals.containsKey(name)) return; // this track is already monitored
        last_intervals.put(name,GenomeLocParser.createGenomeLoc(0,1,1));
        last_rods.put(name,null);
    }

    public static ReferenceOrderedDatum simulate_lookup(String track, GenomeLoc loc, RefMetaDataTracker tracker)  {

        if ( ! last_intervals.containsKey(track)) throw new StingException("Track "+track+" is not monitored by BrokenRODSimulator");

//        if ( loc.getStart() >= 10168704 && loc.getStop() <= 10168728) System.out.println("Request at "+loc);

        ReferenceOrderedDatum oldStyleRod = null; // we will be searching for a record among all the records at a site

        RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData(track,null); // get all records at the site
//        if ( loc.getStart() >= 10168704 && loc.getStop() <= 10168728) {
//            System.out.println("  Rods:" );
//            for ( ReferenceOrderedDatum d : rods ) System.out.println("    "+d.getLocation());
//            System.out.println("  Last ROD is: "+last_intervals.get(track));
//        }

        if ( rods == null || rods.size() == 0 ) return oldStyleRod; // no data, nothing to do

        ReferenceOrderedDatum firstRod = rods.getRecords().get(0);

        // There were quite a few pecularities with the old rod system. First, if there was an "extended" rod
        // (length > 1), and if we landed on it exactly at its start location, we would see that same rod at every
        // reference position until we walk past that rod's stop position. Other rods within the span of that extended rod
        // would be masked (never seen at all). However, if the first time we land inside an extended rod after its start position, we would not
        // see it at all.

        if ( last_intervals.get(track).equals( firstRod.getLocation() ) ) {
            // normally, we would see the first rod spanning current position (can be extended);
            // here we are just making sure that we legitimately "grabbed" this first rod earlier, i.e.
            // landed on its start position
//            if ( loc.getStart() >= 10168704 && loc.getStop() <= 10168728) System.out.println("Returning last");
            return last_rods.get(track);
        }

        // if we are here, the first rod we see at the current location is not the same as the last one we returned
        // in this case we want to skip all extended rods that started before the current position (if any).

        for( ReferenceOrderedDatum d : rods.getRecords() ) {

            if ( d.getLocation().compareTo(loc) < 0 ) continue; // rod starts before current location, old RODIterator would not see it
            oldStyleRod = d;
            break;
        }
        if ( oldStyleRod != null ) {
            last_rods.put(track, oldStyleRod);
            last_intervals.put(track, oldStyleRod.getLocation()); // remember what we just read; note that we would remember an extended rod here only if stepped into it at its start position!
        }
        return oldStyleRod;

    }
}

package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.ReadFeatureAggregator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 */
public class RFWindow {

    private RFAArgumentCollection argumentCollection; // the RF Argument Collection from the RF Walker
    private List<GenomeLoc> windowStarts; // holds the starting positions of the windows
    private Map<String,Boolean> sampleCohortMap; // identifies case samples ("true") and control samples ("false")
    List<ReadFeatureAggregator> aggregators; // read feature aggregators to be used (maintained for fast cloning, and ordering)
    // feature ---> sample ---> count
    List<Map<String,List<ReadFeatureAggregator>>> aggregatorWindows; // holds the map between samples and features in the current windows
    private GenomeLocParser parser; // the active parser, for creating GenomeLocs for empty windows
    private GenomeLoc previousLoc; // the previous genome loc within the user interval given to the RFWindow

    /**
     * Defines a new window which maps samples to read feature aggregators
     * @return - a new element for aggregatorWindows
     */
    private Map<String,List<ReadFeatureAggregator>> newWindow() {
        Map<String,List<ReadFeatureAggregator>> win = new HashMap<String,List<ReadFeatureAggregator>>(sampleCohortMap.size());
        for ( String s : sampleCohortMap.keySet() ) {
            win.put(s,getAggregators());
        }
        // todo -- generalize me
        win.put("case",getAggregators());
        win.put("control",getAggregators());

        return win;
    }

    /**
     * Generates a list of new aggregators to collect data
     * @return list of ReadFeatureAggregators to be the value of a new window
     */
    private List<ReadFeatureAggregator> getAggregators() {
        ArrayList<ReadFeatureAggregator> newEmptyAgs = new ArrayList<ReadFeatureAggregator>(aggregators.size());
        try {
            for ( ReadFeatureAggregator ag : aggregators ) {
                newEmptyAgs.add(ag.getClass().getConstructor(RFAArgumentCollection.class).newInstance(argumentCollection));
            }
        } catch (Exception e) {
            throw new StingException("Error instantiating read feature aggregator",e);
        }

        return newEmptyAgs;
    }

    /**
     * A constructor for the RFWindow object
     * @param collection - the argument collection, from the walker
     * @param cohortMap - the map between sample IDs and whether or not it is a case sample, from the walker
     * @param parser - the Genome Loc Parser, from the walker
     */
    public RFWindow(List<ReadFeatureAggregator> aggregators, RFAArgumentCollection collection, Map<String,Boolean> cohortMap, GenomeLocParser parser) {
        this.argumentCollection = collection;
        this.sampleCohortMap = cohortMap;
        this.parser = parser;
        this.aggregators = aggregators;
    }

    /**
     * Instantiate the tiled windows of sample -> List<aggregator> maps, filling in empty windows up to the one which contains the
     * provided genomeloc if necessary.
     * @param loc - the location to fill up to, usually the starting position of a read or an interval
     * @param currentUserInterval - the current user-provided interval being processed by the walker
     */
    @Requires({"! currentUserInterval.isBefore(loc)"})
    public void instantiate(GenomeLoc loc, GenomeLoc currentUserInterval) {
        aggregatorWindows = new ArrayList<Map<String,List<ReadFeatureAggregator>>>(argumentCollection.windowSize/argumentCollection.windowJump);
        windowStarts = new ArrayList<GenomeLoc>(argumentCollection.windowSize/argumentCollection.windowJump);
        //System.out.printf("Calling fill at %s\t%s%n",loc,currentUserInterval);
        this.fill(loc, currentUserInterval);
    }

    /**
     * Fills the tiled window lists with empty windows up to the one which includes the locus
     * todo -- this can take a lot of memory for large intervals instantiated far into them, e.g.
     *     |---------------------------------------------------------------------| interval
     *                                                                  . <=== locus
     * todo -- perhaps a reduced representation for empty windows that were not created but already have expired?
     * @param loc - the location to fill up to, usually the starting position of a read or an interval
     * @param currentUserInterval - the current user-provided interval being processed by the walker
     */
    @Requires({"! currentUserInterval.isBefore(loc)"})
    public void fill(GenomeLoc loc, GenomeLoc currentUserInterval) {
        // case -1: window could be empty
        if ( windowStarts == null ) {
            instantiate(loc,currentUserInterval);
        }

        // case 0: if the windows are empty add in the beginning of the interval
        if ( windowStarts.size() == 0 ) {
            windowStarts.add(currentUserInterval.getStartLocation());
            aggregatorWindows.add(newWindow());
        }

        // case 1: loc is before or within windowJump bases of the current user interval; we need only instantiate the first window
        if ( loc.isBefore(currentUserInterval) || loc.distance(currentUserInterval) <= argumentCollection.windowJump ) {
            // do nothing at all
        } else {
        // case 2: loc is somewhere in the middle or at the end of current user interval, need to fill in windows up until then
            GenomeLoc nextLoc = windowStarts.get(windowStarts.size()-1);
            while ( loc.distance(nextLoc) > argumentCollection.windowJump ) {
                //System.out.printf("Filling with nextloc %s%n",nextLoc);
                nextLoc =  shiftLocByJump(windowStarts.get(windowStarts.size()-1),argumentCollection.windowJump);
                windowStarts.add(nextLoc);
                aggregatorWindows.add(newWindow());
            }
        }

    }

    /**
     * Shifts a location from chrZ:X to chrZ:X+jump
     * @param loc - location to shift
     * @param jump - amount to shift by
     * @return loc with start position shifted by jump
     */
    @Requires("loc.getStart()+jump < parser.getContigInfo(loc.getContig()).getSequenceLength()") // e.g. that the new loc is not off the end of the contig
    private GenomeLoc shiftLocByJump(GenomeLoc loc, int jump) {
        return parser.createGenomeLoc(loc.getContig(),loc.getStart()+jump);
    }

    /**
     * Fills in missing windows between the previously-seen one and the ending interval, then expires all windows, returning
     * them as complete, resetting previousLocation to null so that the next read re-instantiates internal window data
     * @param userIntervalEnding - the user interval which is ending
     * @return those currently active windows, plus empty interpolating windows between the last active window and the end of the interval
     */
    public List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> flush(GenomeLoc userIntervalEnding) {
        // jump in locations -- flush the windows
        List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> complete = new ArrayList<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>>(aggregators.size());
        // fill in any uncovered windows
        GenomeLoc iStop = userIntervalEnding.getStopLocation();
        //System.out.printf("Calling fill from within flush at %s%n",userIntervalEnding);
        fill(iStop,userIntervalEnding);
        // now expire all windows, terminating them either at the proper window size, or at the endpoint of the interval if it comes sooner
        while ( windowStarts.size() > 0 ) {
            Map<String,List<ReadFeatureAggregator>> cMap = aggregatorWindows.remove(0);
            GenomeLoc wsLoc = windowStarts.remove(0);

            GenomeLoc cLoc;
            if ( wsLoc.distance(iStop) > argumentCollection.windowSize ) {
                cLoc = parser.createGenomeLoc(wsLoc.getContig(),wsLoc.getStart(),wsLoc.getStart()+argumentCollection.windowSize);
            } else {
                cLoc = wsLoc.endpointSpan(iStop);
            }

            complete.add(new Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>(cLoc,cMap));
        }

        previousLoc = null; // will re-instantiate data upon next loc

        return complete;
    }

    /**
     * Workhorse method:
     *    - determines if new windows need be made
     *    - determines if old windows need be tested & trashed
     *    - determines which windows need to be updated & updates them
     *
     * @param loc - starting alignment position of the read
     * @param sample - the sample ID from whom the read was sequenced
     * @return - those feature window(s) that need be tested (and then baleeted)
     */
    public List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> inc(SAMRecord record, GenomeLoc loc, String sample, GenomeLoc userInterval) {
        List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> complete = new ArrayList<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>>(aggregators.size());
        if ( previousLoc == null ) {
            // first time, gotta instantiate stuff
            instantiate(loc,userInterval);
            windowInc(sample,record,aggregatorWindows.get(0));
        } else if ( loc.distance(previousLoc) == 0 ) {
            // best and most common case: just update the living windows
            for ( Map<String,List<ReadFeatureAggregator>> win : aggregatorWindows ) {
                windowInc(sample,record,win);
            }
        } else {
            // another standard case: we've gone to some further base in the same interval
            while ( loc.distance(windowStarts.get(windowStarts.size()-1)) > argumentCollection.windowJump ) {
                // careful, don't use the location itself, but add in windows every winJump bases until this condition is not met
                //System.out.printf("Adding within inc at %s\t%s%n",loc,userInterval);
                windowStarts.add(shiftLocByJump(windowStarts.get(windowStarts.size()-1),argumentCollection.windowJump));
                aggregatorWindows.add(newWindow());
            }

            while ( windowStarts.size() > 0 && loc.distance(windowStarts.get(0)) > argumentCollection.windowSize ) {
                Map<String,List<ReadFeatureAggregator>> cMap = aggregatorWindows.remove(0);
                GenomeLoc wsLoc = windowStarts.remove(0);
                GenomeLoc iStop = userInterval.getStopLocation();
                GenomeLoc cLoc;
                if ( wsLoc.distance(iStop) > argumentCollection.windowSize ) {
                    cLoc = parser.createGenomeLoc(wsLoc.getContig(),wsLoc.getStart(),wsLoc.getStart()+argumentCollection.windowSize);
                } else {
                    cLoc = wsLoc.endpointSpan(iStop);
                }
                complete.add(new Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>(cLoc,cMap));
            }

            for ( Map<String,List<ReadFeatureAggregator>> win : aggregatorWindows ) {
                windowInc(sample,record,win);
            }

        }

        previousLoc = loc;

        return complete;
    }

    /**
     * Incorporate new features into the window
     * @param sample - id of sample from which the features come
     * @param record - the read
     * @param window - the particular window to be updated
     */
    private void windowInc(String sample, SAMRecord record, Map<String,List<ReadFeatureAggregator>> window) {
        if ( sample == null || record == null ) { return; }

        for (ReadFeatureAggregator aggregator : window.get(sample) ) {
            aggregator.aggregate(record);
        }

        for ( ReadFeatureAggregator aggregator : window.get( sampleCohortMap.get(sample) ? "case" : "control") ) {
            aggregator.aggregate(record);
        }
    }

}

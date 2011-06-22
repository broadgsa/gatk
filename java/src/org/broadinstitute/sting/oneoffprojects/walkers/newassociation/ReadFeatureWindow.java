package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features.ReadFeatureAggregator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.*;

/**
 * Workhorse class of read feature association: maintains the active windows, culls the inactive ones and passes
 * them back for testing. Windows consist of aggregators, and list indeces link windows to their starting loc which
 * determines their active/inactive status.
 */
public class ReadFeatureWindow {
    RFAArgumentCollection args;
    private List<GenomeLoc> winStarts;
    private GenomeLoc previousLoc;
    private int winSize;
    private int winJump;
    private GenomeLocParser parser;

    private Map<String,Boolean> sampleCohortMap;
    List<ReadFeatureAggregator> aggregators;
    // feature ---> sample ---> count
    List<Map<String,List<ReadFeatureAggregator>>> aggregatorWindows;

    public ReadFeatureWindow(List<ReadFeatureAggregator> aggregators, Map<String,Boolean> cohortMap, RFAArgumentCollection args, GenomeLocParser parser) {
        winSize = args.windowSize;
        winJump = args.windowJump;
        sampleCohortMap = cohortMap;
        this.args = args;
        this.aggregators =aggregators;
        this.parser = parser;
    }

    public void instantiate(GenomeLoc loc) {
        aggregatorWindows = new ArrayList<Map<String,List<ReadFeatureAggregator>>>(winSize/winJump);
        winStarts = new ArrayList<GenomeLoc>(winSize/winJump);
        winStarts.add(loc);
        aggregatorWindows.add(newWindow());
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
    public List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> inc(SAMRecord record, GenomeLoc loc, String sample) {
        List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> complete = new ArrayList<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>>(aggregators.size());
        if ( previousLoc == null ) {
            // first time, gotta instantiate stuff
            instantiate(loc);
            windowInc(sample,record,aggregatorWindows.get(0));
        } else if ( loc.distance(previousLoc) == 0 ) {
            // best and most common case: just update the living windows
            for ( Map<String,List<ReadFeatureAggregator>> win : aggregatorWindows ) {
                windowInc(sample,record,win);
            }
        } else {
            // another standard case: we've gone to some further base in the same interval
            while ( loc.distance(winStarts.get(winStarts.size()-1)) > winJump ) {
                // careful, don't use the location itself, but add in windows every winJump bases until this condition is not met
                winStarts.add(shiftLocByJump(winStarts.get(winStarts.size()-1),winJump));
                aggregatorWindows.add(newWindow());
            }

            while ( winStarts.size() > 0 && loc.distance(winStarts.get(0)) > winSize ) {
                Map<String,List<ReadFeatureAggregator>> cMap = aggregatorWindows.remove(0);
                GenomeLoc cLoc = winStarts.remove(0).endpointSpan(previousLoc);
                complete.add(new Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>(cLoc,cMap));
            }

            for ( Map<String,List<ReadFeatureAggregator>> win : aggregatorWindows ) {
                windowInc(sample,record,win);
            }

        }

        previousLoc = loc;

        return complete;
    }

    public List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> flush() {
        // jump in locations -- flush the windows
        List<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>> complete = new ArrayList<Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>>(aggregators.size());
        while ( winStarts.size() > 0 ) {
            Map<String,List<ReadFeatureAggregator>> cMap = aggregatorWindows.remove(0);
            GenomeLoc cLoc = winStarts.remove(0).endpointSpan(previousLoc);
            complete.add(new Pair<GenomeLoc,Map<String,List<ReadFeatureAggregator>>>(cLoc,cMap));
        }

        previousLoc = null; // will re-instantiate data upon next loc

        return complete;
    }

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

    private List<ReadFeatureAggregator> getAggregators() {
        ArrayList<ReadFeatureAggregator> newEmptyAgs = new ArrayList<ReadFeatureAggregator>(aggregators.size());
        try {
            for ( ReadFeatureAggregator ag : aggregators ) {
                newEmptyAgs.add(ag.getClass().getConstructor(RFAArgumentCollection.class).newInstance(args));
            }
        } catch (Exception e) {
            throw new StingException("Error instantiating read feature aggregator",e);
        }

        return newEmptyAgs;
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

    private GenomeLoc shiftLocByJump(GenomeLoc loc, int jump) {
        return parser.createGenomeLoc(loc.getContig(),loc.getStart()+jump);
    }
}

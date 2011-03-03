package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/2/11
 * Time: 11:58 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AssociationContext<X,Y> {

    protected List<Map<Sample,Y>> window;

    public AssociationContext() {
        window = new ArrayList<Map<Sample,Y>>(getWindowSize());
    }

    // specifies size of window
    public abstract int getWindowSize();

    // specifies how many bases to wait until next test
    public abstract int slideByValue();

    // specifies whether to use previously seen reads
    public abstract boolean usePreviouslySeenReads();

    // specifies the map from a sample's pileup to the data we want to test
    public abstract Object map(ReadBackedPileup rbp);

    // specifies how to take the per-sample data and reduce them into testable pairs
    public abstract Map<?,X> reduce(List<Map<Sample,Y>> win);

    // do we filter the current location (e.g. omit from window)
    public boolean filter(MapExtender m) { return true; }

    // a basic initialization of the context (give the walker for access to object?)
    public void init(RegionalAssociationWalker walker) { }

    public Map<Sample,Object> mapLocus(MapExtender extender) {
        Map<Sample,ReadBackedPileup> pileups;
        if ( ! usePreviouslySeenReads() ) {
            pileups = extender.getReadFilteredPileup();
        } else {
            pileups = extender.getFullPileup();
        }
        Map<Sample,Object> maps = new HashMap<Sample,Object>(pileups.size());
        for ( Map.Entry<Sample,ReadBackedPileup> samPileup : pileups.entrySet() ) {
            maps.put(samPileup.getKey(),map(samPileup.getValue()));
        }

        return maps;
    }


    public void addData(Map<Sample,Y> sampleData) {
        window.add(sampleData);
    }


    public boolean isFull() {
        return window.size() >= getWindowSize();
    }

    public void slide() {
        ArrayList<Map<Sample,Y>> newWindow = new ArrayList<Map<Sample,Y>>((window.subList(slideByValue(),window.size())));
        newWindow.ensureCapacity(getWindowSize());
        window = newWindow;
    }
}

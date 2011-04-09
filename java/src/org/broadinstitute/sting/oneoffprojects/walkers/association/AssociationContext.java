package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/2/11
 * Time: 11:58 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AssociationContext<X,Y> {

    protected List<Map<Sample,Y>> window;
    private int size;
    private int slide;

    public AssociationContext() {
    }

    public AssociationContext(final RegionalAssociationWalker parent) {
        this.init(parent);
    }

    // specifies whether to use previously seen reads
    public abstract boolean usePreviouslySeenReads();

    // specifies the map from a sample's pileup to the data we want to test
    public abstract Object map(ReadBackedPileup rbp);

    // specifies how to take the per-sample data and reduce them into testable pairs
    public abstract Map<?,X> reduce(List<Map<Sample,Y>> win);

    // do we filter the current location (e.g. omit from window)
    public boolean filter(MapExtender m) { return true; }

    // a basic initialization of the context (give the walker for access to object?)
    public void init(RegionalAssociationWalker walker) {
        size = walker.windowSize;
        slide = walker.slideBy;
        window = new ArrayList<Map<Sample,Y>>(size);
    }

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
        return window.size() >= size;
    }

    public void slide() {
        ArrayList<Map<Sample,Y>> newWindow = new ArrayList<Map<Sample,Y>>((window.subList(slide,window.size())));
        newWindow.ensureCapacity(size);
        window = newWindow;
    }

    public int getWindowSize() {
        return window.size();
    }
}

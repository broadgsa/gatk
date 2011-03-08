package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;

/**
 * @Author chartl
 * @Date 2011-02-23
 * A general context for windowed association
 */
public class RegionalAssociationHandler {
    private MapExtender maps;
    // todo -- the correct way to do this is via the PluginManager (a la VariantEval) but this is a QND implementation
    private Set<AssociationContext> associations;

    public RegionalAssociationHandler(Set<AssociationContext> contexts) {
        maps = new MapExtender();
        associations = contexts;
    }

    public void updateExtender(MapHolder mapHolder) {
        maps.set(mapHolder);
    }

    /**
     * Once the instances are collected; run map-reduce
     * @map: MapExtender --> A child of AssociationContextAtom
     *  @implementation: Indirect construction via newinstance
     * @reduce: (List<AssociationContext>,AssociationContext) --> List<AssociationContextAtom>
     *  @implementation: just append
     */
    public void runMapReduce() {
        for ( AssociationContext w : associations ) {
            if ( w.filter(maps) ) {
                w.addData(w.mapLocus(maps));
            }
        }
    }

    /**
     * Switches what formatting to use based on wiggle or standard
     * @param wiggleFormat - use wiggle format (just Q) or standard (S: P: Q:)
     * @return - test results in proper format
     */
    public Map<AssociationContext,String> runTests(boolean wiggleFormat) {
        if ( wiggleFormat ) {
            return runWiggleTests();
        } else {
            return runTests();
        }
    }

    public Map<AssociationContext,String> runWiggleTests() {
        // todo -- maybe the tdf should be the whole window rather than just the most recent loc?
        Map<AssociationContext,String> testResults = new HashMap<AssociationContext,String>(associations.size());
        for ( AssociationContext w : associations ) {
            if ( w.isFull() ) {
                testResults.put(w,String.format("%d",AssociationTestRunner.getQValue(w)));
                w.slide();
            }
        }
        return testResults;
    }

    /**
     * For those AssociationContexts with full windows:
     *  1) Run their associated test(s)
     *  2) Slide the windows
     */
    public Map<AssociationContext,String> runTests() {
        // todo -- maybe the tdf should be the whole window rather than just the most recent loc?
        Map<AssociationContext,String> testResults = new HashMap<AssociationContext,String>(associations.size());
        for ( AssociationContext w : associations ) {
            if ( w.isFull() ) {
                testResults.put(w,String.format("%s\t%d\t%d\t%s",maps.getReferenceContext().getLocus().getContig(),
                        maps.getReferenceContext().getLocus().getStart(),maps.getReferenceContext().getLocus().getStart()+1,AssociationTestRunner.runTests(w)));
                w.slide();
            }
        }
        return testResults;
    }

    public GenomeLoc getLocation() {
        return maps.getReferenceContext().getLocus();
    }

}
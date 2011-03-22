package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

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
     * @param bedGraphFormat - use bedgraph format (s p q) or standard (S: s P: p Q: q)
     * @return - test results in proper format
     */
    public Map<AssociationContext,String> runTests(boolean bedGraphFormat) {
        // todo -- maybe the tdf should be the whole window rather than just the most recent loc?
        Map<AssociationContext,String> testResults = new HashMap<AssociationContext,String>(associations.size());
        for ( AssociationContext w : associations ) {
            if ( w.isFull() ) {
                String outVal;
                if ( bedGraphFormat ) {
                    Pair<Double,Pair<Double,Integer>> vals = AssociationTestRunner.getTestValues(w);
                    outVal = String.format("%.2f\t%.2e\t%d",vals.first,vals.second.first,vals.second.second);
                } else {
                    outVal = AssociationTestRunner.runTests(w);
                }
                testResults.put(w,String.format("%s\t%d\t%d\t%s",maps.getReferenceContext().getLocus().getContig(),
                        maps.getReferenceContext().getLocus().getStart()-w.getWindowSize()-1,maps.getReferenceContext().getLocus().getStart()+1, outVal));
                w.slide();
            }
        }
        return testResults;
    }

    public GenomeLoc getLocation() {
        return maps.getReferenceContext().getLocus();
    }

}
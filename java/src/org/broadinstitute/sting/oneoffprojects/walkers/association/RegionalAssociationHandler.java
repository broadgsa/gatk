package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
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
    private boolean bedGraphFormat;

    public RegionalAssociationHandler(Set<AssociationContext> contexts, Set<Sample> samples, boolean bedGraph) {
        maps = new MapExtender(samples);
        associations = contexts;
        bedGraphFormat = bedGraph;
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
    public Map<AssociationContext,String> runMapReduce() {
        Map<AssociationContext,String> testResults = new HashMap<AssociationContext,String>();
        if ( maps.getPreviousRef() != null && maps.getPreviousRef().getLocus().compareTo(getLocation()) != -1 ) {
            testResults.putAll(flush());
        }
        for ( AssociationContext w : associations ) {
            if ( w.filter(maps) ) {
                w.addData(w.mapLocus(maps));
            }

            if ( w.isFull() ) {
                testResults.put(w,runTest(w));
                w.slide();
            }
        }

        return testResults;
    }

    /**
     * Switches what formatting to use based on wiggle or standard
     * implicit param bedGraphFormat - use bedgraph format (s p q) or standard (S: s P: p Q: q)
     * @return - test results in proper format
     */
    public String runTest(AssociationContext context) {
        // todo -- maybe the tdf should be the whole window rather than just the most recent loc?
        String outVal;
        if ( bedGraphFormat ) {
            Pair<Double,Pair<Double,Integer>> statVals = AssociationTestRunner.getTestValues(context);
            Pair<Double,Double> simpleDichotVals = AssociationTestRunner.getDichotomizedValues(context);
            outVal = String.format("%.2f\t%.2e\t%d\t%.2f\t%.2f",statVals.first,statVals.second.first,statVals.second.second,
                    simpleDichotVals.first,simpleDichotVals.second);
        } else {
            outVal = AssociationTestRunner.runTests(context);
            Pair<Double,Double> simpleDichotVals = AssociationTestRunner.getDichotomizedValues(context);
            outVal += String.format("\tD: %.2f\tLogD: %.2f",simpleDichotVals.first,simpleDichotVals.second);
        }
        return String.format("%s\t%d\t%d\t%s",maps.getReferenceContext().getLocus().getContig(),
                maps.getReferenceContext().getLocus().getStart()-context.getWindowSize()-1,maps.getReferenceContext().getLocus().getStart()+1, outVal);

    }

    public GenomeLoc getLocation() {
        return maps.getReferenceContext().getLocus();
    }

    /**
     * Flushes context on a jump between intervals. Can not return null.
     * @return on-flush tests
     */
    public Map<AssociationContext,String> flush() {
        Map<AssociationContext,String> flushTests = new HashMap<AssociationContext,String>();
        for ( AssociationContext ac : associations ) {
            if ( ac.canFlush() ) {
                if ( ac.testOnFlush() ) {
                    flushTests.put(ac,runTest(ac));
                }
                ac.flush();
            }
        }

        return flushTests;
    }

    public Set<AssociationContext> getAssociations() {
        return associations;
    }

    public MapExtender getExtender() {
        return maps;
    }
}
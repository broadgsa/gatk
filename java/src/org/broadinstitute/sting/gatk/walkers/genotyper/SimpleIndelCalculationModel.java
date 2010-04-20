package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;

import java.util.*;

public class SimpleIndelCalculationModel extends GenotypeCalculationModel {

    // the previous normal event context
    private Map<String, StratifiedAlignmentContext> cachedContext;

    protected SimpleIndelCalculationModel() {}

    private int testSkipCount = 5;

    public VariantCallContext callLocus(RefMetaDataTracker tracker, char ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors) {
        cachedContext = contexts;
        return null;
    }

    public VariantCallContext callExtendedLocus(RefMetaDataTracker tracker, char[] ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts) {

        System.out.println("\nReached " + loc + " through an extended event");
        for (Map.Entry<String,StratifiedAlignmentContext> e : contexts.entrySet()) {
            System.out.println("Set "+e.getKey());
            System.out.println("  Context: "+e.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size()+ " reads");
            System.out.println("  Cached context: "+
                      cachedContext.get(e.getKey()).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size()+ " reads");
            System.out.println("  First read in cached context: "+
                      cachedContext.get(e.getKey()).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().getReads().get(0).getReadName());
            ReadBackedExtendedEventPileup p = e.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup();
            if ( p== null ) System.out.println("EXTENDED PILEUP IS NULL");
            System.out.println("  Event(s): " + e.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup().getEventStringsWithCounts(ref));
        }
        if ( testSkipCount==0 ) System.exit(1);
        testSkipCount--;
        // TODO -- implement me

        //VariantContext vc = new MutableVariantContext("UG_indel_call", loc, alleles, genotypes, phredScaledConfidence/10.0, null, attributes);
        //return new VariantCallContext(vc, phredScaledConfidence >= CONFIDENCE_THRESHOLD);

        return null;
    }
}
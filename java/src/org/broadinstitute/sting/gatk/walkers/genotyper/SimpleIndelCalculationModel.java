package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.*;

import java.util.*;

public class SimpleIndelCalculationModel extends GenotypeCalculationModel {

    private int MIN_COVERAGE = 6;
    private double MIN_FRACTION = 0.3;
    private double MIN_CONSENSUS_FRACTION = 0.7 ;
    // the previous normal event context
//    private Map<String, StratifiedAlignmentContext> cachedContext;

    protected SimpleIndelCalculationModel() {}

    private int totalIndels = 0;
    private int totalCoverage = 0;
    private int bestIndelCount = 0;
    private String bestEvent = null;


    public VariantCallContext callLocus(RefMetaDataTracker tracker, byte ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors) {
//        cachedContext = contexts;
        return null;
    }

    public VariantCallContext callExtendedLocus(RefMetaDataTracker tracker, byte[] ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts) {

        totalIndels = 0;
        totalCoverage = 0;
        bestIndelCount = 0;
        bestEvent = null;
 /*
        System.out.println("\nReached " + loc + " through an extended event");
        for (Map.Entry<String,StratifiedAlignmentContext> e : contexts.entrySet()) {
            System.out.println("Set "+e.getKey());
            System.out.println("  Context: "+e.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size()+ " reads");
            ReadBackedExtendedEventPileup p = e.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup();
            if ( p== null ) System.out.println("EXTENDED PILEUP IS NULL");
            System.out.println("  Event(s): " + e.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup().getEventStringsWithCounts(ref));
        }
*/
        initializeAlleles(ref, contexts);

        VariantCallContext vcc = new VariantCallContext(false);

        if ( totalIndels == 0 ) return vcc; // this can happen if indel-containing reads get filtered out by the engine

        if ( totalCoverage < MIN_COVERAGE ) return vcc;

        if ( ((double)bestIndelCount)/totalCoverage < MIN_FRACTION ) return vcc;

        if ( ((double)bestIndelCount)/totalIndels < MIN_CONSENSUS_FRACTION ) return vcc;

        List<Allele> alleles = new ArrayList<Allele>(2);

        if ( bestEvent.charAt(0) == '+') {
            alleles.add( Allele.create(Allele.NULL_ALLELE_STRING,true) );
            alleles.add( Allele.create(bestEvent.substring(1), false ));
        } else {
            if ( bestEvent.charAt(0) == '-' ) {
                alleles.add( Allele.create(Allele.NULL_ALLELE_STRING,false) );
                alleles.add( Allele.create(bestEvent.substring(1), true ));
                loc = GenomeLocParser.setStop(loc, loc.getStop() + bestEvent.length()-1);
            } else
                throw new GATKException("Internal error (probably a bug): event does not conform to expected format: "+ bestEvent);
        }

        VariantContext vc = new VariantContext("UG_Indel_call", loc.getContig(), loc.getStart(), loc.getStop(), alleles, new HashMap<String, Genotype>() /* genotypes */,
                -1.0 /* log error */, null /* filters */, null /* attributes */);

        vcc = new VariantCallContext(vc,true);
        vcc.setRefBase(ref[0]);
/*
        if ( totalIndels > 0 ) {
            System.out.println("Calling: "+bestEvent+" ["+bestIndelCount+"/"+totalIndels+"/"+totalCoverage+"] at "+loc);
        } else {
            System.out.println("NO EVENT");
        }
*/
        //VariantContext vc = new MutableVariantContext("UG_indel_call", loc, alleles, genotypes, phredScaledConfidence/10.0, null, attributes);
        //return new VariantCallContext(vc, phredScaledConfidence >= CONFIDENCE_THRESHOLD);

        return vcc;
    }

    protected void initializeAlleles(byte[] ref, Map<String, StratifiedAlignmentContext> contexts) {


        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            AlignmentContext context = sample.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

            totalCoverage += context.size();

            // calculate the sum of quality scores for each base
            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup();

            List<Pair<String,Integer>> all_events = pileup.getEventStringsWithCounts(ref);
            for ( Pair<String,Integer> p : all_events ) {
                if ( p.second > bestIndelCount ) {
                    bestIndelCount = p.second;
                    bestEvent = p.first;
                }
                totalIndels += p.second;
            }
        }
    }
}
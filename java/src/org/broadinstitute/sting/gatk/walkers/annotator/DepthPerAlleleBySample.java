package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFFormatHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.pileup.*;

import java.util.*;


public class DepthPerAlleleBySample implements GenotypeAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, StratifiedAlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
        if ( g == null || !g.isCalled() )
            return null;

        if ( !vc.isBiallelic() )
            return null;

        if ( vc.isSNP() )
            return annotateSNP(stratifiedContext, vc);
        if ( vc.isIndel() )
            return annotateIndel(stratifiedContext, vc);

        return null;
    }

    private Map<String,Object> annotateSNP(StratifiedAlignmentContext stratifiedContext, VariantContext vc) {

        HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
        alleleCounts.put(vc.getReference().getBases()[0], 0);
        alleleCounts.put(vc.getAlternateAllele(0).getBases()[0], 0);

        ReadBackedPileup pileup = stratifiedContext.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
        for ( PileupElement p : pileup ) {
            if ( alleleCounts.containsKey(p.getBase()) )
                alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+1);
        }

        // we need to add counts in the correct order
        Integer[] counts = new Integer[2];
        counts[0] = alleleCounts.get(vc.getReference().getBases()[0]);
        counts[1] = alleleCounts.get(vc.getAlternateAllele(0).getBases()[0]);

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), counts);
        return map;
    }

    private Map<String,Object> annotateIndel(StratifiedAlignmentContext stratifiedContext, VariantContext vc) {
        ReadBackedExtendedEventPileup pileup = stratifiedContext.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup();
        if ( pileup == null )
            return null;

        Integer[] counts = new Integer[2];

        // TODO -- fix me
        /*
        for ( ExtendedEventPileupElement e : pileup.toExtendedIterable() ) {
            if ( countsBySize.keySet().contains(e.getEventLength()) ) { // if proper length
                if ( e.isDeletion() && vc.isDeletion() || e.isInsertion() && vc.isInsertion() ) {
                    countsBySize.put(e.getEventLength(),countsBySize.get(e.getEventLength())+1);
                }
            }
        }
        */

        Map<String,Object> map = new HashMap<String,Object>();
        map.put(getKeyNames().get(0), counts);
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("AD"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), 2, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles; currently only works for bi-allelic sites")); }
}
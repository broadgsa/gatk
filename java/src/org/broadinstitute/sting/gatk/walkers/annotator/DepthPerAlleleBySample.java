package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFFormatHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.pileup.*;

import java.util.*;


public class DepthPerAlleleBySample implements GenotypeAnnotation, ExperimentalAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, StratifiedAlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
        if ( g == null || !g.isCalled() )
            return null;

        if ( vc.isSNP() ) {
            return annotateSNP(stratifiedContext, vc);
        } else if ( vc.isIndel() ) {
            return annotateIndel(stratifiedContext, vc);
        } else {
            return null;
        }
    }

    private Map<String,Object> annotateSNP(StratifiedAlignmentContext stratifiedContext, VariantContext vc) {

        Set<Allele> altAlleles = vc.getAlternateAlleles();
        if ( altAlleles.size() == 0 )
            return null;

        HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
        for ( Allele allele : altAlleles )
            alleleCounts.put(allele.getBases()[0], 0);

        ReadBackedPileup pileup = stratifiedContext.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
        for (PileupElement p : pileup ) {
            if ( alleleCounts.containsKey(p.getBase()) )
                alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+1);
        }

        StringBuffer sb = new StringBuffer();
        // we need to add counts in the correct order
        for ( Allele allele : vc.getAlleles() ) {
            if ( alleleCounts.containsKey(allele.getBases()[0]) ) {
                if ( sb.length() > 0 )
                    sb.append(",");
                sb.append(String.format("%d", alleleCounts.get(allele.getBases()[0])));
            }
        }

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), sb.toString());
        return map;
    }

    private Map<String,Object> annotateIndel(StratifiedAlignmentContext stratifiedContext, VariantContext vc) {
        ReadBackedExtendedEventPileup pileup = stratifiedContext.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getExtendedEventPileup();
        if ( pileup == null ) {
            return null;
        }
        // get identities and lengths for indel events
        // TODO -- Insertions and deletions represented at the same locus
        HashMap<Integer,Integer> countsBySize = new HashMap<Integer,Integer>();
        for ( Allele al : vc.getAlternateAlleles() ) {
            countsBySize.put(al.length(),0);
        }

        for ( ExtendedEventPileupElement e : pileup.toExtendedIterable() ) {
            if ( countsBySize.keySet().contains(e.getEventLength()) ) { // if proper length
                if ( e.isDeletion() && vc.isDeletion() || e.isInsertion() && vc.isInsertion() ) {
                    countsBySize.put(e.getEventLength(),countsBySize.get(e.getEventLength())+1);
                }
            }
        }

        StringBuffer sb = new StringBuffer();
        char type = vc.isDeletion() ? 'D' : 'I';

        for ( int len : countsBySize.keySet() ) {
            if ( sb.length() > 0 ) {
                sb.append(',');
            }
            sb.append(String.format("%d%s%d",len,type,countsBySize.get(len)));
        }

        Map<String,Object> map = new HashMap<String,Object>();
        map.put(getKeyNames().get(0),sb.toString());
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("AD"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), -1, VCFHeaderLineType.Integer, "Depth in genotypes for each ALT allele, in the same order as listed")); }
}
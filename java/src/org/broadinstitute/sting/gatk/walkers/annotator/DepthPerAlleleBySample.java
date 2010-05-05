package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.*;


public class DepthPerAlleleBySample implements GenotypeAnnotation, ExperimentalAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, StratifiedAlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
        // for now, we don't support indels
        if ( g == null || !g.isCalled() || vc.getType() != VariantContext.Type.SNP )
            return null;

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
        map.put(getKeyName(), sb.toString());
        return map;
    }

    public String getKeyName() { return "AD"; }

    public VCFFormatHeaderLine getDescription() { return new VCFFormatHeaderLine(getKeyName(), 1, VCFFormatHeaderLine.FORMAT_TYPE.Integer, "Depth in genotypes for each ALT allele, in the same order as listed"); }
}
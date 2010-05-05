package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Map;
import java.util.HashMap;


public class HomopolymerRun implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        if ( !vc.isBiallelic() || !vc.isSNP() )
            return null;

        int run = computeHomopolymerRun(vc.getAlternateAllele(0).toString().charAt(0), ref);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%d", run));
        return map;
    }

    public String getKeyName() { return "HRun"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("HRun", 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Largest Contiguous Homopolymer Run of Variant Allele In Either Direction"); }

    public boolean useZeroQualityReads() { return false; }

    private static int computeHomopolymerRun(char altAllele, ReferenceContext ref) {

        // TODO -- this needs to be computed in a more accurate manner
        // We currently look only at direct runs of the alternate allele adjacent to this position

        char[] bases = ref.getBases();
        GenomeLoc window = ref.getWindow();
        GenomeLoc locus = ref.getLocus();

        int refBasePos = (int)(locus.getStart() - window.getStart());

        int leftRun = 0;
        for ( int i = refBasePos - 1; i >= 0; i--) {
            if ( Character.toUpperCase(bases[i]) != altAllele )
                break;
            leftRun++;
        }

        int rightRun = 0;
        for ( int i = refBasePos + 1; i < bases.length; i++) {
            if ( Character.toUpperCase(bases[i]) != altAllele )
                break;
            rightRun++;
        }

        return Math.max(leftRun, rightRun);
     }
}
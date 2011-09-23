package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Largest contiguous homopolymer run of the variant allele in either direction on the reference.
 */
public class HomopolymerRun extends InfoFieldAnnotation implements StandardAnnotation {

    private boolean ANNOTATE_INDELS = true;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {

        if ( !vc.isBiallelic() )
            return null;

        int run;
        if ( vc.isSNP() ) {
            run = computeHomopolymerRun(vc.getAlternateAllele(0).getBases()[0], ref);
        } else if ( vc.isIndel() && ANNOTATE_INDELS ) {
            run = computeIndelHomopolymerRun(vc,ref);
        } else {
            return null;
        }
        
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", run));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("HRun"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("HRun", 1, VCFHeaderLineType.Integer, "Largest Contiguous Homopolymer Run of Variant Allele In Either Direction")); }

    public boolean useZeroQualityReads() { return false; }

    private static int computeHomopolymerRun(byte altAllele, ReferenceContext ref) {

        // TODO -- this needs to be computed in a more accurate manner
        // We currently look only at direct runs of the alternate allele adjacent to this position

        byte[] bases = ref.getBases();
        GenomeLoc window = ref.getWindow();
        GenomeLoc locus = ref.getLocus();

        int refBasePos = (int)(locus.getStart() - window.getStart());

        int leftRun = 0;
        for ( int i = refBasePos - 1; i >= 0; i--) {
            if ( bases[i] != altAllele )
                break;
            leftRun++;
        }

        int rightRun = 0;
        for ( int i = refBasePos + 1; i < bases.length; i++) {
            if ( bases[i] != altAllele )
                break;
            rightRun++;
        }

        return Math.max(leftRun, rightRun);
     }

    private static int computeIndelHomopolymerRun(VariantContext vc, ReferenceContext ref) {
        byte[] bases = ref.getBases();
        GenomeLoc locus = ref.getLocus();
        GenomeLoc window = ref.getWindow();
        int refBasePos = (int) (locus.getStart() - window.getStart())+1;
        if ( vc.isSimpleDeletion() ) {
            // check that deleted bases are the same
            byte dBase = bases[refBasePos];
            for ( int i = 0; i < vc.getReference().length(); i ++ ) {
                if ( bases[refBasePos+i] != dBase ) {
                    return 0;
                }
            }

            return computeHomopolymerRun(dBase, ref);
        } else {
            // check that inserted bases are the same
            byte insBase = vc.getAlternateAllele(0).getBases()[0];
            for ( byte b : vc.getAlternateAllele(0).getBases() ) {
                if ( insBase != (char) b ) {
                    return 0;
                }
            }

            return computeHomopolymerRun(insBase,ref);
        }
    }
}
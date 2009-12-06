package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.Variation;


public class HomopolymerRun extends StandardVariantAnnotation {

    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation) {

        if ( !variation.isBiallelic() || !variation.isSNP() )
            return null;

        int run = computeHomopolymerRun(variation.getAlternativeBaseForSNP(), ref);
        return String.format("%d", run);
    }

    public String getKeyName() { return "HRun"; }

    public String getDescription() { return "HRun,1,Integer,\"Largest Contiguous Homopolymer Run of Variant Allele In Either Direction\""; }

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
}
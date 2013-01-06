package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 1/3/13
 * Time: 11:36 AM
 * To change this template use File | Settings | File Templates.
 */
public class AverageAltAlleleLength extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation, ExperimentalAnnotation {

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Average Allele Length"));
    }

    public List<String> getKeyNames() { return Arrays.asList("AAL"); }

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        if ( !vc.hasLog10PError() )
            return null;

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        Map<String, Object> map = new HashMap<String, Object>();

        double length = getMeanAltAlleleLength(vc);
        map.put(getKeyNames().get(0),String.format("%.2f",length));
        return map;
    }

    public static double getMeanAltAlleleLength(VariantContext vc) {
        double averageLength = 1.0;
        if ( ! vc.isSNP() && ! vc.isSymbolic() ) {
            // adjust for the event length
            int averageLengthNum = 0;
            int averageLengthDenom = 0;
            int refLength = vc.getReference().length();
            for ( Allele a : vc.getAlternateAlleles() ) {
                int numAllele = vc.getCalledChrCount(a);
                int alleleSize;
                if ( a.length() == refLength ) {
                    // SNP or MNP
                    byte[] a_bases = a.getBases();
                    byte[] ref_bases = vc.getReference().getBases();
                    int n_mismatch = 0;
                    for ( int idx = 0; idx < a_bases.length; idx++ ) {
                        if ( a_bases[idx] != ref_bases[idx] )
                            n_mismatch++;
                    }
                    alleleSize = n_mismatch;
                }
                else if ( a.isSymbolic() ) {
                    alleleSize = 1;
                } else {
                    alleleSize = Math.abs(refLength-a.length());
                }
                averageLengthNum += alleleSize*numAllele;
                averageLengthDenom += numAllele;
            }
            averageLength = ( (double) averageLengthNum )/averageLengthDenom;
        }

        return averageLength;
    }
}

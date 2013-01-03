package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Variant confidence (from the QUAL field) / unfiltered depth.
 *
 * Low scores are indicative of false positive calls and artifacts.  Note that QualByDepth requires sequencing
 * reads associated with the samples with polymorphic genotypes.
 */
public class QualByDepth extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

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

        int depth = 0;

        for ( final Genotype genotype : genotypes ) {

            // we care only about variant calls with likelihoods
            if ( !genotype.isHet() && !genotype.isHomVar() )
                continue;

            if (stratifiedContexts!= null) {
                AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());
                if ( context == null )
                    continue;
                depth += context.getBasePileup().depthOfCoverage();

            }
            else if (perReadAlleleLikelihoodMap != null) {
                PerReadAlleleLikelihoodMap perReadAlleleLikelihoods = perReadAlleleLikelihoodMap.get(genotype.getSampleName());
                if (perReadAlleleLikelihoods == null || perReadAlleleLikelihoods.isEmpty())
                    continue;

                depth += perReadAlleleLikelihoods.getNumberOfStoredElements();
            }
        }

        if ( depth == 0 )
            return null;

        double QD = -10.0 * vc.getLog10PError() / (double)depth;
        Map<String, Object> map = new HashMap<String, Object>();

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
            double averageLength = ( (double) averageLengthNum )/averageLengthDenom;
            QD /= averageLength;
            map.put(getKeyNames().get(1),String.format("%.2f",averageLength));
        }

        map.put(getKeyNames().get(0), String.format("%.2f", QD));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("QD","AAL"); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth"),
                             new VCFInfoHeaderLine(getKeyNames().get(1), 1, VCFHeaderLineType.Float, "Average Allele Length"));
    }


}

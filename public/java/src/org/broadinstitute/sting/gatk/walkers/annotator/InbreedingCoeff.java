package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;


/**
 * Likelihood-based (using PL field) test for the inbreeding among samples.
 *
 * A continuous generalization of the Hardy-Weinberg test for disequilibrium that works
 * well with limited coverage per sample.  See the 1000 Genomes Phase I release for
 * more information.  Note that the Inbreeding Coefficient will not be calculated for files
 * with fewer than a minimum (generally 10) number of samples.
 */
public class InbreedingCoeff extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    private static final int MIN_SAMPLES = 10;
    private Set<String> founderIds;

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        if(founderIds == null && walker != null)
            founderIds =  ((Walker)walker).getSampleDB().getFounderIds();
        return calculateIC(vc);
    }

    private Map<String, Object> calculateIC(final VariantContext vc) {
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
        if ( genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant())
            return null;

        int idxAA = 0, idxAB = 1, idxBB = 2;

        if (!vc.isBiallelic()) {
            // for non-bliallelic case, do test with most common alt allele.
            // Get then corresponding indeces in GL vectors to retrieve GL of AA,AB and BB.
            int[] idxVector = vc.getGLIndecesOfAlternateAllele(vc.getAltAlleleWithHighestAlleleCount());
            idxAA = idxVector[0];
            idxAB = idxVector[1];
            idxBB = idxVector[2];
        }

        double refCount = 0.0;
        double hetCount = 0.0;
        double homCount = 0.0;
        int N = 0; // number of samples that have likelihoods
        for ( final Genotype g : genotypes ) {
            if ( g.isNoCall() || !g.hasLikelihoods() )
                continue;

            if (g.getPloidy() != 2) // only work for diploid samples
                continue;

            N++;
            final double[] normalizedLikelihoods = MathUtils.normalizeFromLog10( g.getLikelihoods().getAsVector() );
            refCount += normalizedLikelihoods[idxAA];
            hetCount += normalizedLikelihoods[idxAB];
            homCount += normalizedLikelihoods[idxBB];
        }

        if( N < MIN_SAMPLES ) {
            return null;
        }

        final double p = ( 2.0 * refCount + hetCount ) / ( 2.0 * (refCount + hetCount + homCount) ); // expected reference allele frequency
        final double q = 1.0 - p; // expected alternative allele frequency
        final double F = 1.0 - ( hetCount / ( 2.0 * p * q * (double)N ) ); // inbreeding coefficient

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.4f", F));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("InbreedingCoeff"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("InbreedingCoeff", 1, VCFHeaderLineType.Float, "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation")); }
}
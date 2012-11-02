package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Trio;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.RodRequiringAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Given a variant context, uses the genotype likelihoods to assess the likelihood of the site being a mendelian violation
 * versus the likelihood of the site transmitting according to mendelian rules. This assumes that the organism is
 * diploid. When multiple trios are present, the annotation is simply the maximum of the likelihood ratios, rather than
 * the strict 1-Prod(1-p_i) calculation, as this can scale poorly for uncertain sites and many trios.
 */

public class MVLikelihoodRatio extends InfoFieldAnnotation implements ExperimentalAnnotation, RodRequiringAnnotation {

    private MendelianViolation mendelianViolation = null;
    public static final String MVLR_KEY = "MVLR";
    private Set<Trio> trios;

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( mendelianViolation == null ) {
            trios = ((Walker) walker).getSampleDB().getTrios();
            if ( trios.size() > 0 ) {
                mendelianViolation = new MendelianViolation(((VariantAnnotator)walker).minGenotypeQualityP );
            }
            else {
                throw new UserException("Mendelian violation annotation can only be used from the Variant Annotator, and must be provided a valid PED file (-ped) from the command line.");
            }
        }

        Map<String,Object> attributeMap = new HashMap<String,Object>(1);
        //double pNoMV = 1.0;
        double maxMVLR = Double.MIN_VALUE;
        for ( Trio trio : trios ) {
            if ( contextHasTrioLikelihoods(vc,trio) ) {
                Double likR = mendelianViolation.violationLikelihoodRatio(vc,trio.getMaternalID(),trio.getPaternalID(),trio.getChildID());
                maxMVLR = likR > maxMVLR ? likR : maxMVLR;
                //pNoMV *= (1.0-Math.pow(10.0,likR)/(1+Math.pow(10.0,likR)));
            }
        }

        //double pSomeMV = 1.0-pNoMV;
        //toRet.put("MVLR",Math.log10(pSomeMV)-Math.log10(1.0-pSomeMV));
        if ( Double.compare(maxMVLR,Double.MIN_VALUE) != 0 )
            attributeMap.put(MVLR_KEY,maxMVLR);
        return attributeMap;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList(MVLR_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(MVLR_KEY, 1, VCFHeaderLineType.Float, "Mendelian violation likelihood ratio: L[MV] - L[No MV]")); }


    private boolean contextHasTrioLikelihoods(VariantContext context, Trio trio) {
        for ( String sample : Arrays.asList(trio.getMaternalID(),trio.getPaternalID(),trio.getChildID()) ) {
            if ( ! context.hasGenotype(sample) )
                return false;
            if ( ! context.getGenotype(sample).hasLikelihoods() )
                return false;
        }

        return true;
    }

}

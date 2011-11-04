package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFilterHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/14/11
 * Time: 12:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class MVLikelihoodRatio extends InfoFieldAnnotation implements ExperimentalAnnotation {

    private MendelianViolation mendelianViolation = null;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( mendelianViolation == null ) {
            if ( walker instanceof VariantAnnotator ) {
                mendelianViolation = new MendelianViolation(((VariantAnnotator)walker).familyStr, ((VariantAnnotator)walker).minGenotypeQualityP );
            }
            else {
                throw new UserException("Mendelian violation annotation can only be used from the Variant Annotator");
            }
        }

        Map<String,Object> toRet = new HashMap<String,Object>(1);
        boolean hasAppropriateGenotypes = vc.hasGenotype(mendelianViolation.getSampleChild()) &&
                vc.hasGenotype(mendelianViolation.getSampleDad()) &&
                vc.hasGenotype(mendelianViolation.getSampleMom());
        if ( hasAppropriateGenotypes )
            toRet.put("MVLR",mendelianViolation.violationLikelihoodRatio(vc));

        return toRet;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList("MVLR"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MVLR", 1, VCFHeaderLineType.Float, "Mendelian violation likelihood ratio: L[MV] - L[No MV]")); }
}

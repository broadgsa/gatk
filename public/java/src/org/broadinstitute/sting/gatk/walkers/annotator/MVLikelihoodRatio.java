package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.samples.SampleDB;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.RodRequiringAnnotation;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/14/11
 * Time: 12:24 PM
 */

public class MVLikelihoodRatio extends InfoFieldAnnotation implements ExperimentalAnnotation, RodRequiringAnnotation {

    private MendelianViolation mendelianViolation = null;
    private String motherId;
    private String fatherId;
    private String childId;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( mendelianViolation == null ) {
            if (checkAndSetSamples(((Walker) walker).getSampleDB())) {
                mendelianViolation = new MendelianViolation(((VariantAnnotator)walker).minGenotypeQualityP );
            }
            else {
                throw new UserException("Mendelian violation annotation can only be used from the Variant Annotator, and must be provided a valid PED file (-ped) from the command line containing only 1 trio.");
            }
        }

        Map<String,Object> toRet = new HashMap<String,Object>(1);
        boolean hasAppropriateGenotypes = vc.hasGenotype(motherId) && vc.getGenotype(motherId).hasLikelihoods() &&
                vc.hasGenotype(fatherId) && vc.getGenotype(fatherId).hasLikelihoods() &&
                vc.hasGenotype(childId) && vc.getGenotype(childId).hasLikelihoods();
        if ( hasAppropriateGenotypes )
            toRet.put("MVLR",mendelianViolation.violationLikelihoodRatio(vc,motherId,fatherId,childId));

        return toRet;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList("MVLR"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MVLR", 1, VCFHeaderLineType.Float, "Mendelian violation likelihood ratio: L[MV] - L[No MV]")); }

    private boolean checkAndSetSamples(SampleDB db){
        Set<String> families = db.getFamilyIDs();
        if(families.size() != 1)
            return false;

        Set<Sample> family = db.getFamily(families.iterator().next());
        if(family.size() != 3)
            return false;

        Iterator<Sample> sampleIter = family.iterator();
        Sample sample;
        for(sample = sampleIter.next();sampleIter.hasNext();sample=sampleIter.next()){
            if(sample.getParents().size()==2){
                motherId = sample.getMaternalID();
                fatherId = sample.getPaternalID();
                childId = sample.getID();
                return true;
            }
        }
        return false;
    }

}

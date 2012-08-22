package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.samples.SampleDB;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.RodRequiringAnnotation;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
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
    private Set<Trio> trios;
    private class Trio {
        String motherId;
        String fatherId;
        String childId;
    }

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( mendelianViolation == null ) {
            if (checkAndSetSamples(((Walker) walker).getSampleDB())) {
                mendelianViolation = new MendelianViolation(((VariantAnnotator)walker).minGenotypeQualityP );
            }
            else {
                throw new UserException("Mendelian violation annotation can only be used from the Variant Annotator, and must be provided a valid PED file (-ped) from the command line.");
            }
        }

        Map<String,Object> toRet = new HashMap<String,Object>(1);
        //double pNoMV = 1.0;
        double maxMVLR = Double.MIN_VALUE;
        for ( Trio trio : trios ) {
            boolean hasAppropriateGenotypes = vc.hasGenotype(trio.motherId) && vc.getGenotype(trio.motherId).hasLikelihoods() &&
                    vc.hasGenotype(trio.fatherId) && vc.getGenotype(trio.fatherId).hasLikelihoods() &&
                    vc.hasGenotype(trio.childId) && vc.getGenotype(trio.childId).hasLikelihoods();
            if ( hasAppropriateGenotypes ) {
                Double likR = mendelianViolation.violationLikelihoodRatio(vc,trio.motherId,trio.fatherId,trio.childId);
                maxMVLR = likR > maxMVLR ? likR : maxMVLR;
                //pNoMV *= (1.0-Math.pow(10.0,likR)/(1+Math.pow(10.0,likR)));
            }
        }

        //double pSomeMV = 1.0-pNoMV;
        //toRet.put("MVLR",Math.log10(pSomeMV)-Math.log10(1.0-pSomeMV));
        toRet.put("MVLR",maxMVLR);
        return toRet;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList("MVLR"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MVLR", 1, VCFHeaderLineType.Float, "Mendelian violation likelihood ratio: L[MV] - L[No MV]")); }

    private boolean checkAndSetSamples(SampleDB db){
        trios = new HashSet<Trio>();
        Set<String> families = db.getFamilyIDs();
        for ( String familyString : families ) {
            Set<Sample> family = db.getFamily(familyString);
            Iterator<Sample> sampleIterator = family.iterator();
            Sample sample;
            for ( sample = sampleIterator.next(); sampleIterator.hasNext(); sample=sampleIterator.next()) {
                if ( sample.getParents().size() == 2 ) {
                    Trio trio = new Trio();
                    trio.childId = sample.getID();
                    trio.fatherId = sample.getFather().getID();
                    trio.motherId = sample.getMother().getID();
                    trios.add(trio);
                }
            }
        }

        return trios.size() > 0;
    }

}

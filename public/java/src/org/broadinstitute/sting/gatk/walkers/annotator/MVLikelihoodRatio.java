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
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
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
 * the strict 1-∏(1-p_i) calculation, as this can scale poorly for uncertain sites and many trios.
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
            trios = checkAndSetSamples(((Walker) walker).getSampleDB());
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
                Double likR = mendelianViolation.violationLikelihoodRatio(vc,trio.getMaternalID(),trio.getPaternalID(),trio.childId);
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

    // todo - this entire function should be in samples DB
    private Set<Trio> checkAndSetSamples(SampleDB db){
        Set<Trio> trioSet = new HashSet<Trio>();
        for ( String familyString : db.getFamilyIDs() ) {
            Set<Sample> family = db.getFamily(familyString);
            for ( Sample sample : family) {
                if ( sample.getParents().size() == 2 ) {
                    Trio trio = new Trio(sample.getMaternalID(),sample.getPaternalID(),sample.getID());
                    trioSet.add(trio);
                }
            }
        }

        return trioSet;
    }

    private boolean contextHasTrioLikelihoods(VariantContext context, Trio trio) {
        for ( String sample : trio ) {
            if ( ! context.hasGenotype(sample) )
                return false;
            if ( ! context.getGenotype(sample).hasLikelihoods() )
                return false;
        }

        return true;
    }

    // TODO -- this class is too much.
    // TODO -- Why iterable?
    // TODO -- shuoldn't this be in samplesDB() so you can just called samplesDB().getTrios()
    // TODO -- should just have final string IDs, and getters, no setters
    private class Trio implements Iterable<String> {
        private String maternalID;
        private String paternalID;
        private String childId;

        public Trio(String mom, String dad, String child) {
            this.maternalID = mom;
            this.paternalID = dad;
            this.childId = child;
        }

        public String getMaternalID() {
            return this.maternalID;
        }

        public String getPaternalID() {
            return this.paternalID;
        }

        public String getChildId() {
            return this.childId;
        }

        public void setMaternalID(String id) {
            this.maternalID = id;
        }

        public void setPaternalID(String id) {
            this.paternalID = id;
        }

        public void setChildId(String id) {
            this.childId = id;
        }

        public Iterator<String> iterator() {
            return Arrays.asList(maternalID,paternalID,childId).iterator();
        }
    }
}

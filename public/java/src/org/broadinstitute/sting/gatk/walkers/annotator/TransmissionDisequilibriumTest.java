package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 11/14/11
 */

public class TransmissionDisequilibriumTest extends InfoFieldAnnotation implements ExperimentalAnnotation {

    private Set<MendelianViolation> fullMVSet = null;
    private final static int REF = 0;
    private final static int HET = 1;
    private final static int HOM = 2;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( fullMVSet == null ) {
            fullMVSet = new HashSet<MendelianViolation>();

            if ( walker instanceof VariantAnnotator ) {
                final Map<String,Set<Sample>> families = ((VariantAnnotator) walker).getSampleDB().getFamilies();
                for( final Set<Sample> family : families.values() ) {
                    for( final Sample sample : family ) {
                        if( sample.getParents().size() == 2 && family.containsAll(sample.getParents()) ) { // only works with trios for now
                            fullMVSet.add( new MendelianViolation(sample, 0.0) );
                        }
                    }
                }
            } else {
                throw new UserException("Transmission disequilibrium test annotation can only be used from the Variant Annotator and requires a valid ped file be passed in.");
            }
        }

        final Map<String,Object> toRet = new HashMap<String,Object>(1);
        final HashSet<MendelianViolation> mvsToTest = new HashSet<MendelianViolation>();

        for( final MendelianViolation mv : fullMVSet ) {
            final boolean hasAppropriateGenotypes = vc.hasGenotype(mv.getSampleChild()) && vc.getGenotype(mv.getSampleChild()).hasLikelihoods() &&
                    vc.hasGenotype(mv.getSampleDad()) && vc.getGenotype(mv.getSampleDad()).hasLikelihoods() &&
                    vc.hasGenotype(mv.getSampleMom()) && vc.getGenotype(mv.getSampleMom()).hasLikelihoods();
            if ( hasAppropriateGenotypes ) {
                mvsToTest.add(mv);
            }
        }

        toRet.put("TDT", calculateTDT( vc, mvsToTest ));

        return toRet;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList("TDT"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("TDT", 1, VCFHeaderLineType.Float, "Test statistic from Wittkowski transmission disequilibrium test.")); }

    // Following derivation in http://en.wikipedia.org/wiki/Transmission_disequilibrium_test#A_modified_version_of_the_TDT
    private double calculateTDT( final VariantContext vc, final Set<MendelianViolation> mvsToTest ) {

        final double nABGivenABandBB = calculateNChildren(vc, mvsToTest, HET, HET, HOM);
        final double nBBGivenABandBB = calculateNChildren(vc, mvsToTest, HOM, HET, HOM);
        final double nAAGivenABandAB = calculateNChildren(vc, mvsToTest, REF, HET, HET);
        final double nBBGivenABandAB = calculateNChildren(vc, mvsToTest, HOM, HET, HET);
        final double nAAGivenAAandAB = calculateNChildren(vc, mvsToTest, REF, REF, HET);
        final double nABGivenAAandAB = calculateNChildren(vc, mvsToTest, HET, REF, HET);

        final double numer = (nABGivenABandBB - nBBGivenABandBB) + 2.0 * (nAAGivenABandAB - nBBGivenABandAB) + (nAAGivenAAandAB - nABGivenAAandAB);
        final double denom = (nABGivenABandBB + nBBGivenABandBB) + 4.0 * (nAAGivenABandAB + nBBGivenABandAB) + (nAAGivenAAandAB + nABGivenAAandAB);
        return (numer * numer) / denom;
    }

    private double calculateNChildren( final VariantContext vc, final Set<MendelianViolation> mvsToTest, final int childIdx, final int momIdx, final int dadIdx ) {
        final double likelihoodVector[] = new double[mvsToTest.size() * 2];
        int iii = 0;
        for( final MendelianViolation mv : mvsToTest ) {
            final double[] momGL = vc.getGenotype(mv.getSampleMom()).getLikelihoods().getAsVector();
            final double[] dadGL = vc.getGenotype(mv.getSampleDad()).getLikelihoods().getAsVector();
            final double[] childGL = vc.getGenotype(mv.getSampleChild()).getLikelihoods().getAsVector();
            likelihoodVector[iii++] = momGL[momIdx] + dadGL[dadIdx] + childGL[childIdx];
            likelihoodVector[iii++] = momGL[dadIdx] + dadGL[momIdx] + childGL[childIdx];
        }

        return MathUtils.sumLog10(likelihoodVector);
    }
}

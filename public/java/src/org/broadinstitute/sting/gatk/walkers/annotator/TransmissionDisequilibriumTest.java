package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin, lfran
 * Date: 11/14/11
 */

public class TransmissionDisequilibriumTest extends InfoFieldAnnotation implements ExperimentalAnnotation {

    private Set<Sample> trios = null;
    private final static int REF = 0;
    private final static int HET = 1;
    private final static int HOM = 2;
    private final static int MIN_NUM_VALID_TRIOS = 5; // don't calculate this population-level statistic if there are less than X trios with full genotype likelihood information

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( trios == null ) {
            if ( walker instanceof VariantAnnotator ) {
                trios =  ((VariantAnnotator) walker).getSampleDB().getChildrenWithParents();
            } else {
                throw new UserException("Transmission disequilibrium test annotation can only be used from the Variant Annotator and requires a valid ped file be passed in.");
            }
        }

        final Map<String,Object> toRet = new HashMap<String,Object>(1);
        final HashSet<Sample> triosToTest = new HashSet<Sample>();

        for( final Sample child : trios) {
            final boolean hasAppropriateGenotypes = vc.hasGenotype(child.getID()) && vc.getGenotype(child.getID()).hasLikelihoods() &&
                    vc.hasGenotype(child.getPaternalID()) && vc.getGenotype(child.getPaternalID()).hasLikelihoods() &&
                    vc.hasGenotype(child.getMaternalID()) && vc.getGenotype(child.getMaternalID()).hasLikelihoods();
            if ( hasAppropriateGenotypes ) {
                triosToTest.add(child);
            }
        }

        if( triosToTest.size() >= MIN_NUM_VALID_TRIOS ) {
            toRet.put("TDT", calculateTDT( vc, triosToTest ));
        }

        return toRet;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList("TDT"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("TDT", 1, VCFHeaderLineType.Float, "Test statistic from Wittkowski transmission disequilibrium test.")); }

    // Following derivation in http://en.wikipedia.org/wiki/Transmission_disequilibrium_test#A_modified_version_of_the_TDT
    private double calculateTDT( final VariantContext vc, final Set<Sample> triosToTest ) {

        final double nABGivenABandBB = calculateNChildren(vc, triosToTest, HET, HET, HOM) + calculateNChildren(vc, triosToTest, HET, HOM, HET);
        final double nBBGivenABandBB = calculateNChildren(vc, triosToTest, HOM, HET, HOM) + calculateNChildren(vc, triosToTest, HOM, HOM, HET);
        final double nAAGivenABandAB = calculateNChildren(vc, triosToTest, REF, HET, HET);
        final double nBBGivenABandAB = calculateNChildren(vc, triosToTest, HOM, HET, HET);
        final double nAAGivenAAandAB = calculateNChildren(vc, triosToTest, REF, REF, HET) + calculateNChildren(vc, triosToTest, REF, HET, REF);
        final double nABGivenAAandAB = calculateNChildren(vc, triosToTest, HET, REF, HET) + calculateNChildren(vc, triosToTest, HET, HET, REF);

        final double numer = (nABGivenABandBB - nBBGivenABandBB) + 2.0 * (nAAGivenABandAB - nBBGivenABandAB) + (nAAGivenAAandAB - nABGivenAAandAB);
        final double denom = (nABGivenABandBB + nBBGivenABandBB) + 4.0 * (nAAGivenABandAB + nBBGivenABandAB) + (nAAGivenAAandAB + nABGivenAAandAB);
        return (numer * numer) / denom;
    }

    private double calculateNChildren( final VariantContext vc, final Set<Sample> triosToTest, final int childIdx, final int parent1Idx, final int parent2Idx ) {
        final double likelihoodVector[] = new double[triosToTest.size()];
        int iii = 0;
        for( final Sample child : triosToTest ) {
            final double[] momGL = vc.getGenotype(child.getMaternalID()).getLikelihoods().getAsVector();
            final double[] dadGL = vc.getGenotype(child.getPaternalID()).getLikelihoods().getAsVector();
            final double[] childGL = vc.getGenotype(child.getID()).getLikelihoods().getAsVector();
            likelihoodVector[iii++] = momGL[parent1Idx] + dadGL[parent2Idx] + childGL[childIdx];
        }

        return MathUtils.sumLog10(likelihoodVector);
    }
}

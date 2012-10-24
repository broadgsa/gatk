package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.RodRequiringAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin, lfran, ebanks
 * Date: 11/14/11
 */

public class TransmissionDisequilibriumTest extends InfoFieldAnnotation implements ExperimentalAnnotation, RodRequiringAnnotation {

    private Set<Sample> trios = null;
    private final static int MIN_NUM_VALID_TRIOS = 5; // don't calculate this population-level statistic if there are less than X trios with full genotype likelihood information

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( trios == null ) {
            if ( walker instanceof VariantAnnotator ) {
                trios =  ((VariantAnnotator) walker).getSampleDB().getChildrenWithParents();
            } else {
                throw new UserException("Transmission disequilibrium test annotation can only be used from the Variant Annotator and requires a valid ped file be passed in.");
            }
        }

        final Map<String, Object> toRet = new HashMap<String, Object>(1);
        final HashSet<Sample> triosToTest = new HashSet<Sample>();

        for( final Sample child : trios ) {
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

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("TDT", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Test statistic from Wittkowski transmission disequilibrium test.")); }

    // Following derivation in http://en.wikipedia.org/wiki/Transmission_disequilibrium_test#A_modified_version_of_the_TDT
    private List<Double> calculateTDT( final VariantContext vc, final Set<Sample> triosToTest ) {

        List<Double> pairwiseTDTs = new ArrayList<Double>(10);
        final int HomRefIndex = 0;

        // for each pair of alleles, add the likelihoods
        int numAltAlleles = vc.getAlternateAlleles().size();
        for ( int alt = 1; alt <= numAltAlleles; alt++ ) {
            final int HetIndex = alt;
            final int HomVarIndex = determineHomIndex(alt, numAltAlleles+1);

            final double nABGivenABandBB = calculateNChildren(vc, triosToTest, HetIndex, HetIndex, HomVarIndex) + calculateNChildren(vc, triosToTest, HetIndex, HomVarIndex, HetIndex);
            final double nBBGivenABandBB = calculateNChildren(vc, triosToTest, HomVarIndex, HetIndex, HomVarIndex) + calculateNChildren(vc, triosToTest, HomVarIndex, HomVarIndex, HetIndex);
            final double nAAGivenABandAB = calculateNChildren(vc, triosToTest, HomRefIndex, HetIndex, HetIndex);
            final double nBBGivenABandAB = calculateNChildren(vc, triosToTest, HomVarIndex, HetIndex, HetIndex);
            final double nAAGivenAAandAB = calculateNChildren(vc, triosToTest, HomRefIndex, HomRefIndex, HetIndex) + calculateNChildren(vc, triosToTest, HomRefIndex, HetIndex, HomRefIndex);
            final double nABGivenAAandAB = calculateNChildren(vc, triosToTest, HetIndex, HomRefIndex, HetIndex) + calculateNChildren(vc, triosToTest, HetIndex, HetIndex, HomRefIndex);

            final double numer = (nABGivenABandBB - nBBGivenABandBB) + 2.0 * (nAAGivenABandAB - nBBGivenABandAB) + (nAAGivenAAandAB - nABGivenAAandAB);
            final double denom = (nABGivenABandBB + nBBGivenABandBB) + 4.0 * (nAAGivenABandAB + nBBGivenABandAB) + (nAAGivenAAandAB + nABGivenAAandAB);
            pairwiseTDTs.add((numer * numer) / denom);
        }

        return pairwiseTDTs;
    }

    private double calculateNChildren( final VariantContext vc, final Set<Sample> triosToTest, final int childIdx, final int momIdx, final int dadIdx ) {
        final double likelihoodVector[] = new double[triosToTest.size()];
        int iii = 0;
        for( final Sample child : triosToTest ) {
            final double[] momGL = vc.getGenotype(child.getMaternalID()).getLikelihoods().getAsVector();
            final double[] dadGL = vc.getGenotype(child.getPaternalID()).getLikelihoods().getAsVector();
            final double[] childGL = vc.getGenotype(child.getID()).getLikelihoods().getAsVector();
            likelihoodVector[iii++] = momGL[momIdx] + dadGL[dadIdx] + childGL[childIdx];
        }

        return MathUtils.sumLog10(likelihoodVector);
    }
    
    private static int determineHomIndex(final int alleleIndex, int numAlleles) {
        int result = 0;
        for ( int i = 0; i < alleleIndex; i++ )
            result += numAlleles--;
        return result;
    }
}

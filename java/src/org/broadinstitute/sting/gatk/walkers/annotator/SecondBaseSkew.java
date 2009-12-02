package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 16, 2009
 * Time: 11:25:51 AM
 * To change this template use File | Settings | File Templates.
 */
public class SecondBaseSkew implements VariantAnnotation {
    private final static double epsilon = Math.pow(10.0,-12);
    private final static boolean USE_ZERO_QUALITY_READS = false; // todo -- should be false in my opinion MAD 
    private final static String KEY_NAME = "2b_Chi";
    private final static double[] UNIFORM_ON_OFF_RATIO = {1.0/3, 2.0/3};
    private double[] proportionExpectations = UNIFORM_ON_OFF_RATIO;

    public boolean useZeroQualityReads() { return USE_ZERO_QUALITY_READS; }

    public String getKeyName() { return KEY_NAME; }

    public String getDescription() { return KEY_NAME + ",1,Float,\"Chi-square Secondary Base Skew\""; }

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        if ( variation.isSNP() && variation.isBiallelic() ) {
            char snp = variation.getAlternativeBaseForSNP();
//            try {
//                System.out.printf("snp %c, alt is %c%n", snp, getNonref(genotypes, ref.getBase()));
//            } catch (IllegalStateException e) {
//                System.out.printf("%s is not biallelic%n", variation.toString());
//                return null;
//            }

            Pair<Integer,Double> depthProp = getSecondaryPileupNonrefEstimator(ref.getBase(), pileup, snp);
            if ( depthProp == null ) {
                return null;
            } else {
                //System.out.printf("%d / %f%n", depthProp.getFirst(), depthProp.getSecond());
                double p_transformed = transform(depthProp.getSecond(), depthProp.getFirst());
                double expected_transformed = transform(proportionExpectations[0], depthProp.getFirst());
                // System.out.println("p_transformed="+p_transformed+" e_transformed="+expected_transformed+" variantDepth="+depthProp.getFirst());
                // System.out.println("Proportion variant bases with ref 2bb="+depthProp.getSecond()+" Expected="+proportionExpectations[0]);
                double chi_square =  Math.signum(depthProp.getSecond() - proportionExpectations[0])*Math.min(Math.pow(p_transformed - expected_transformed, 2), Double.MAX_VALUE);
                return new Pair<String,String>(KEY_NAME, String.format("%f", chi_square));
            }
        } else {
            return null;
        }
    }

    private double transform( double proportion, int depth ) {
        proportion = proportion - epsilon;
        return proportion / ( Math.sqrt ( proportion*(1-proportion)/depth ) );
    }

    private Pair<Integer, Double> getSecondaryPileupNonrefEstimator(char ref, ReadBackedPileup p, char snp ) {
        int variantDepth = 0;
        int variantsWithRefSecondBase = 0;

        for (PileupElement pile : p ) {
            byte pbase = pile.getBase();
            byte sbase = pile.getSecondBase();

            if ( BaseUtils.isRegularBase((char)sbase) && BaseUtils.basesAreEqual(pbase, (byte) snp) ) {
                variantDepth++;
                if ( BaseUtils.basesAreEqual(sbase, (byte)ref) ) {
                    variantsWithRefSecondBase++;
                }
            }
        }

        if ( variantDepth > 0 ) {
            //System.out.printf("%d %d %d %d%n", primaryPileup.length, secondaryPileup.length, variantDepth, variantsWithRefSecondBase );
            double biasedProportion = ( 1.0 + variantsWithRefSecondBase )/(1.0 + variantDepth );
            return new Pair<Integer,Double>(variantDepth+1, biasedProportion);
        } else {
            return null;
        }
    }

    private char getNonref(List<Genotype> genotypes, char ref) {
        for ( Genotype g : genotypes ) {
            if ( g.isVariant(ref) ) {
                return g.toVariation(ref).getAlternativeBaseForSNP();
            }
        }

        throw new IllegalStateException("List of genotypes did not contain a variant.");
    }
}

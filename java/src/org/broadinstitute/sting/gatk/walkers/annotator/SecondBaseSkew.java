package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 16, 2009
 * Time: 11:25:51 AM
 * To change this template use File | Settings | File Templates.
 */
public class SecondBaseSkew implements VariantAnnotation{
    private static double epsilon = Math.pow(10.0,-12);
    private static boolean USE_ZERO_QUALITY_READS = true;
    private static String KEY_NAME = "2b_Chi";
    private static double[] UNIFORM_ON_OFF_RATIO = {1.0/3, 2.0/3};
    private double[] proportionExpectations = UNIFORM_ON_OFF_RATIO;


    public boolean useZeroQualityReads() { return USE_ZERO_QUALITY_READS; }

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {
        double chi_square;
        Pair<Integer,Double> depthProp = getSecondaryPileupNonrefEstimator(pileup,genotypes);
            if ( depthProp == null ) {
                return null;
            } else {
                double p_transformed = transform(depthProp.getSecond(), depthProp.getFirst());
                double expected_transformed = transform(proportionExpectations[0], depthProp.getFirst());
                // System.out.println("p_transformed="+p_transformed+" e_transformed="+expected_transformed+" variantDepth="+depthProp.getFirst());
                // System.out.println("Proportion variant bases with ref 2bb="+depthProp.getSecond()+" Expected="+proportionExpectations[0]);
                chi_square =  Math.signum(depthProp.getSecond() - proportionExpectations[0])*Math.min(Math.pow(p_transformed - expected_transformed, 2), Double.MAX_VALUE);
            }

        return new Pair<String,String>(KEY_NAME,String.format("%f", chi_square));
    }

    private double transform( double proportion, int depth ) {
        proportion = proportion - epsilon;
        return proportion / ( Math.sqrt ( proportion*(1-proportion)/depth ) );
    }

    private Pair<Integer, Double> getSecondaryPileupNonrefEstimator(ReadBackedPileup p, List<Genotype> genotypes) {
        char snp;
        try {
            snp = getNonref(genotypes, p.getRef());
        } catch ( IllegalStateException e ) {
            // tri-allelic site
            // System.out.println("Illegal State Exception caught at "+p.getLocation().toString()+" 2bb skew annotation suppressed ("+e.getLocalizedMessage()+")");
            return null;
        }
        int variantDepth = 0;
        int variantsWithRefSecondBase = 0;
        char[] primaryPileup = p.getBases().toCharArray();
        String secondBasePileup = p.getSecondaryBasePileup();

        if ( secondBasePileup == null ) {
            // System.out.println("Warning: Second base pileup is null at "+p.getLocation().toString());
            return null;
        } else {
            char [] secondaryPileup = secondBasePileup.toCharArray();
            //System.out.printf("primary=%d secondary=%d locus=%s%n", primaryPileup.length, secondaryPileup.length, p.getLocation().toString());

            for ( int i = 0; i < primaryPileup.length; i ++ ) {
                if ( BaseUtils.basesAreEqual((byte) primaryPileup[i], (byte) snp) ) {
                    variantDepth++;
                    if ( BaseUtils.basesAreEqual((byte) secondaryPileup[i], (byte) p.getRef()) ) {
                        variantsWithRefSecondBase++;
                    }
                }
            }


            double biasedProportion = ( 1.0 + variantsWithRefSecondBase )/(1.0 + variantDepth );

            return new Pair<Integer,Double>(variantDepth+1, biasedProportion);
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

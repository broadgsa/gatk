package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.apache.log4j.Logger;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 21, 2009
 * Time: 6:08:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class PrimaryBaseSecondaryBaseSymmetry implements VariantAnnotation{

    private static boolean USE_ZERO_MAPQ_READS = true;
    private static boolean USE_ZERO_QUALITY_READS = false;
    private static String KEY_NAME = "1b2b_symmetry";
    private static double epsilon = Math.pow(10.0,-12);
    Logger logger = Logger.getLogger(PrimaryBaseSecondaryBaseSymmetry.class);

    private boolean useConservativeVariance = true;

    public void conservativeVarianceUsage( boolean b ) {
        useConservativeVariance = b;
    }

    public boolean useZeroMappingQualityReads() { return USE_ZERO_MAPQ_READS; }

    public boolean useZeroQualityReads() { return USE_ZERO_QUALITY_READS; }

    public Pair<String,String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {
        Pair<Integer,Double> refSecondBasePair = getProportionOfReferenceSecondBasesThatSupportAlt(ref, pileup, genotypes);
        Pair<Integer,Double> nonrefPrimaryBasePair = getProportionOfPrimaryNonrefBasesThatSupportAlt(ref, pileup, genotypes);
        if ( refSecondBasePair == null || nonrefPrimaryBasePair ==  null ) {
            return null;
        } else {
            logger.info("2b="+refSecondBasePair.second+" 1b="+nonrefPrimaryBasePair.second);
            //double primary_secondary_stat = 1.0/Math.pow( transform(refSecondBasePair) - transform(nonrefPrimaryBasePair), 2);
            double primary_secondary_stat = refSecondBasePair.second - nonrefPrimaryBasePair.second;
            String annotation = String.format("%f", primary_secondary_stat);
            logger.info("Second-base symmetry: annotating with "+annotation);
            return new Pair<String,String>(KEY_NAME, annotation);
        }
    }

    private double transform( double proportion, int depth ) {
        proportion = proportion - epsilon;
        if ( useConservativeVariance ) {
            return proportion / ( Math.sqrt ( 0.5*(1-0.5) / Math.sqrt(depth) ) );
        } else {
            return proportion / ( Math.sqrt ( proportion*(1-proportion)/depth ) );
        }
    }

    private double transform( Pair<Integer, Double> depth_prop ) {
        return transform( depth_prop.getSecond(), depth_prop.getFirst() );
    }

    private Pair<Integer,Double> getProportionOfReferenceSecondBasesThatSupportAlt( ReferenceContext ref, ReadBackedPileup p, List<Genotype> genotypes) {
        char snp;
        try {
            snp = getNonref(genotypes, ref.getBase());
        } catch ( IllegalStateException e) {
            logger.info("Caught: IllegalStateException -- "+e.getLocalizedMessage());
            return null;
        }

        String secondaryPileupStr = p.getSecondaryBasePileup();
        if ( secondaryPileupStr == null )
            return null;

        char[] secondaryPileup = secondaryPileupStr.toCharArray();
        int depth = p.size();

        int support = 0;
        for ( char c : secondaryPileup ) {
            if ( Character.toUpperCase(c) == Character.toUpperCase(snp) ) {
                support ++;
            }
        }

        double as_prop = ( ( double ) support ) / depth;

        return new Pair<Integer,Double> ( depth, as_prop );
    }

    private Pair<Integer,Double> getProportionOfPrimaryNonrefBasesThatSupportAlt( ReferenceContext ref, ReadBackedPileup p, List<Genotype> genotypes ) {
        char snp;
        try {
            snp = getNonref(genotypes, ref.getBase());
        } catch ( IllegalStateException e ) {
            return null;
        }

        int [] baseCounts = p.getBasePileupAsCounts();
        int support = -1;
        int depth = 0;
        for ( char c : BaseUtils.BASES ) {
            // ignore ref
            if ( Character.toUpperCase(c) == Character.toUpperCase(ref.getBase()) ) {
            } else {
            // catch our snp
                if ( Character.toUpperCase(c) == Character.toUpperCase(snp) ) {
                    support = baseCounts[BaseUtils.simpleBaseToBaseIndex(c)];
                    depth = depth + baseCounts[BaseUtils.simpleBaseToBaseIndex(c)];
                } else {
                    depth = depth + baseCounts[BaseUtils.simpleBaseToBaseIndex(c)];
                }
            }
        }

        if ( depth == 0 || support < 0 ) {
            return null;
        }

        double as_prop = ( ( double ) support) / depth;

        return new Pair<Integer,Double> ( depth, as_prop );
    }

    private char getNonref(List<Genotype> genotypes, char ref) {
        //logger.info(genotypes.size());
        for ( Genotype g : genotypes ) {
            //logger.info("Genotype: "+g.getBases()+" Ref from genotype: "+g.getReference()+" Ref from method: "+ref);
            if ( g.isVariant(ref) ) {
                return g.toVariation(ref).getAlternativeBaseForSNP();
            }
        }
        throw new IllegalStateException("List of genotypes did not contain a variant.");
    }

}

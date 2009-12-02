package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
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
    //
    // Where are the integration tests for this piece of code?
    //
    private static boolean USE_ZERO_MAPQ_READS = false;
    private static String KEY_NAME = "1b2b_symmetry";
    Logger logger = Logger.getLogger(PrimaryBaseSecondaryBaseSymmetry.class);

    public boolean useZeroQualityReads() { return USE_ZERO_MAPQ_READS; }

    public Pair<String,String> annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        // todo -- this code doesn't work, should't be called
        if ( true )
            return null;
        else {
            if ( variation.isSNP() && variation.isBiallelic() ) {
                byte snp = (byte)variation.getAlternativeBaseForSNP();
                Pair<Integer,Double> refSecondBasePair = getProportionOfReferenceSecondBasesThatSupportAlt(ref, pileup, snp);
                Pair<Integer,Double> nonrefPrimaryBasePair = getProportionOfPrimaryNonrefBasesThatSupportAlt(ref, pileup, (char)snp);
                if ( refSecondBasePair == null || nonrefPrimaryBasePair ==  null ) {
                    return null;
                } else {
                    //System.out.printf("refSecondBasePair = %s, nonrefPrimaryBasePair = %s%n", refSecondBasePair, nonrefPrimaryBasePair);
                    double primary_secondary_stat = refSecondBasePair.second - nonrefPrimaryBasePair.second;
                    String annotation = String.format("%f", primary_secondary_stat);
                    return new Pair<String,String>(KEY_NAME, annotation);
                }
            } else {
                return null;
            }
        }
    }

    private Pair<Integer,Double> getProportionOfReferenceSecondBasesThatSupportAlt( ReferenceContext ref, ReadBackedPileup p, byte snp ) {
        int depth = 0;
        int support = 0;
        byte refBase = (byte)ref.getBase();

        for (PileupElement pile : p ) {
            byte c = pile.getSecondBase();

            if ( BaseUtils.isRegularBase(c) && BaseUtils.basesAreEqual(pile.getBase(), refBase)) { // stops indels et al
                depth++;
                support += BaseUtils.basesAreEqual(c, snp) ? 1 : 0;
            }
        }

        if ( depth > 0 ) {
            double as_prop = ( ( double ) support ) / depth;
            return new Pair<Integer,Double> ( depth, as_prop );
        } else {
            return null;
        }
    }

    private Pair<Integer,Double> getProportionOfPrimaryNonrefBasesThatSupportAlt( ReferenceContext ref, ReadBackedPileup p, char snp ) {
        // todo -- Why is it looping?
        int [] baseCounts = p.getBaseCounts();
        int support = -1;
        int depth = 0;
        for ( char c : BaseUtils.BASES ) {
            // ignore ref
            if ( Character.toUpperCase(c) != Character.toUpperCase(ref.getBase()) ) {
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
}

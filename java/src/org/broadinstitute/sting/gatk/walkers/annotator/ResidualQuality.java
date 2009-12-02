package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 23, 2009
 * Time: 1:39:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class ResidualQuality implements VariantAnnotation{
    private static double EPSILON = Math.pow(10,-12);
    public static String KEY_NAME = "ResidualQuality";

    public boolean useZeroQualityReads() { return true; } // for robustness

    public Pair<String,String> annotate( ReferenceContext ref, ReadBackedPileup p, Variation variation, List<Genotype> genotypes) {

        Character snp = getSNPChar(ref, genotypes);
        if ( snp == null ) {
            return null;
        }

        Double logResidQual = getLogResidualQuality(p,ref.getBase(),snp);

        if ( logResidQual == null ) {
            return null;
        }

        return new Pair<String,String>(KEY_NAME, String.format("%f", logResidQual ));
    }

    public String getKeyName() { return KEY_NAME; }

    public String getDescription() { return KEY_NAME + ",1,Float,\"Log-scaled Residual Error\""; }

    private Double getLogResidualQuality( ReadBackedPileup p, char ref, char snp ) {
        byte[] pbp = p.getBases();
        byte[] quals = p.getQuals();
        if ( pbp == null || quals == null ) {
            return null;
        }

        int nSNP = 0;
        int nRef = 0;
        int pileupSize = pbp.length;
        int sumOfQuals = 0;
        for ( int i = 0; i < pileupSize; i ++ ) {
            if ( BaseUtils.basesAreEqual( pbp[i], (byte) ref ) ) {
                // ref site
                nRef ++;
            } else if ( BaseUtils.basesAreEqual ( pbp[i], ( byte ) snp ) ) {
                // snp site
                nSNP++;
            } else {
                // non-ref non-snp site, increase quality
                sumOfQuals += quals[i];
            }
        }

        // want to return sum of qualities times the proportion of non-SNP bases they account for (includes ref bases)
        // taking the log gives    log(SQ) + log(non-ref non-snp bases) - log(non-snp bases)

        return Math.log(sumOfQuals + EPSILON) + Math.log(pileupSize-nSNP-nRef+EPSILON) - Math.log(pileupSize-nSNP + EPSILON);
    }

    private Character getSNPChar( ReferenceContext ref, List<Genotype> genotypes ) {
        try {
            return getNonref( genotypes, ref.getBase() );
        } catch ( IllegalStateException e ) {
            return null;
        }
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

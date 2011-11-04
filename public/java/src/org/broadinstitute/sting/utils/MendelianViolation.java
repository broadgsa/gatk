package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Collection;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * User: carneiro
 * Date: 3/9/11
 * Time: 12:38 PM
 */
public class MendelianViolation {
    String sampleMom;
    String sampleDad;
    String sampleChild;

    List allelesMom;
    List allelesDad;
    List allelesChild;

    double minGenotypeQuality;

    static final int[] mvOffsets = new int[] { 1,2,5,6,8,11,15,18,20,21,24,25 };
    static final int[] nonMVOffsets = new int[]{ 0,3,4,7,9,10,12,13,14,16,17,19,22,23,26 };

    private static Pattern FAMILY_PATTERN = Pattern.compile("(.*)\\+(.*)=(.*)");

    public String getSampleMom() {
        return sampleMom;
    }
    public String getSampleDad() {
        return sampleDad;
    }
    public String getSampleChild() {
        return sampleChild;
    }
    public double getMinGenotypeQuality() {
        return minGenotypeQuality;
    }

    /**
     *
     * @param sampleMomP - sample name of mom
     * @param sampleDadP - sample name of dad
     * @param sampleChildP - sample name of child
     */
    public MendelianViolation (String sampleMomP, String sampleDadP, String sampleChildP) {
        sampleMom = sampleMomP;
        sampleDad = sampleDadP;
        sampleChild = sampleChildP;
    }

    /**
     *
     * @param family - the sample names string "mom+dad=child"
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     */
    public MendelianViolation(String family, double minGenotypeQualityP) {
        minGenotypeQuality = minGenotypeQualityP;

        Matcher m = FAMILY_PATTERN.matcher(family);
        if (m.matches()) {
            sampleMom = m.group(1);
            sampleDad = m.group(2);
            sampleChild = m.group(3);
        }
        else
            throw new IllegalArgumentException("Malformatted family structure string: " + family + " required format is mom+dad=child");
    }

    /**
     * An alternative to the more general constructor if you want to get the Sample information from the engine yourself.
     * @param sample - the sample object extracted from the sample metadata YAML file given to the engine.
     * @param minGenotypeQualityP - the minimum phred scaled genotype quality score necessary to asses mendelian violation
     */
    public MendelianViolation(Sample sample, double minGenotypeQualityP) {
        sampleMom = sample.getMother().getID();
        sampleDad = sample.getFather().getID();
        sampleChild = sample.getID();
        minGenotypeQuality = minGenotypeQualityP;
    }

    /**
     * This method prepares the object to evaluate for violation. Typically you won't call it directly, a call to
     * isViolation(vc) will take care of this. But if you want to know whether your site was a valid comparison site
     * before evaluating it for mendelian violation, you can call setAlleles and then isViolation().
     * @param vc - the variant context to extract the genotypes and alleles for mom, dad and child.
     * @return false if couldn't find the genotypes or context has empty alleles. True otherwise.
     */
    public boolean setAlleles (VariantContext vc)
    {
        Genotype gMom = vc.getGenotypes(sampleMom).get(sampleMom);
        Genotype gDad = vc.getGenotypes(sampleDad).get(sampleDad);
        Genotype gChild = vc.getGenotypes(sampleChild).get(sampleChild);

        if (gMom == null || gDad == null || gChild == null)
            throw new IllegalArgumentException(String.format("Variant %s:%d didn't contain genotypes for all family members: mom=%s dad=%s child=%s", vc.getChr(), vc.getStart(), sampleMom, sampleDad, sampleChild));

        if (gMom.isNoCall() || gDad.isNoCall() || gChild.isNoCall() ||
            gMom.getPhredScaledQual()   < minGenotypeQuality ||
            gDad.getPhredScaledQual()   < minGenotypeQuality ||
            gChild.getPhredScaledQual() < minGenotypeQuality ) {

            return false;
        }

        allelesMom = gMom.getAlleles();
        allelesDad = gDad.getAlleles();
        allelesChild = gChild.getAlleles();
        return !allelesMom.isEmpty() && !allelesDad.isEmpty() && !allelesChild.isEmpty();
    }


    /**
     *
     * @param vc the variant context to extract the genotypes and alleles for mom, dad and child.
     * @return False if we can't determine (lack of information), or it's not a violation. True if it is a violation.
     *
     */
    public boolean isViolation(VariantContext vc)
    {
        return setAlleles(vc) && isViolation();
    }

    /**
     * @return whether or not there is a mendelian violation at the site.
     */
    public boolean isViolation() {
        if (allelesMom.contains(allelesChild.get(0)) && allelesDad.contains(allelesChild.get(1)) ||
            allelesMom.contains(allelesChild.get(1)) && allelesDad.contains(allelesChild.get(0)))
            return false;
        return true;
    }

    /**
     * @return the likelihood ratio for a mendelian violation
     */
    public double violationLikelihoodRatio(VariantContext vc) {
        double[] logLikAssignments = new double[27];
        // the matrix to set up is
        // MOM   DAD    CHILD
        //                    |-  AA
        //   AA     AA    |    AB
        //                    |-   BB
        //                    |- AA
        //  AA     AB     |   AB
        //                    |- BB
        // etc. The leaves are counted as 0-11 for MVs and 0-14 for non-MVs
        double[] momGL = vc.getGenotype(sampleMom).getLikelihoods().getAsVector();
        double[] dadGL = vc.getGenotype(sampleDad).getLikelihoods().getAsVector();
        double[] childGL = vc.getGenotype(sampleChild).getLikelihoods().getAsVector();
        int offset = 0;
        for ( int oMom = 0; oMom < 3; oMom++ ) {
            for ( int oDad = 0; oDad < 3; oDad++ ) {
                for ( int oChild = 0; oChild < 3; oChild ++ ) {
                    logLikAssignments[offset++] = momGL[oMom] + dadGL[oDad] + childGL[oChild];
                }
            }
        }
        double[] mvLiks = new double[12];
        double[] nonMVLiks = new double[15];
        for ( int i = 0; i < 12; i ++ ) {
            mvLiks[i] = logLikAssignments[mvOffsets[i]];
        }

        for ( int i = 0; i < 15; i++) {
            nonMVLiks[i] = logLikAssignments[nonMVOffsets[i]];
        }

        return MathUtils.log10sumLog10(mvLiks) - MathUtils.log10sumLog10(nonMVLiks);
    }

}

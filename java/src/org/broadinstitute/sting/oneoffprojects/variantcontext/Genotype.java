package org.broadinstitute.sting.oneoffprojects.variantcontext;

import java.util.*;

/**
 * @author ebanks
 *         <p/>
 *         Class Genotype
 *         <p/>
 *         This class emcompasses all the basic information about a genotype
 */
public class Genotype extends AttributedObject {
    private List<Allele> alleles;

    private double negLog10PError;

    private String sample;

    public Genotype(List<Allele> alleles, String sample, double negLog10PError) {
        this.alleles = new ArrayList<Allele>(alleles);
        this.sample = sample;
        this.negLog10PError = negLog10PError;
    }

    /**
     * @return the alleles for this genotype
     */
    public List<Allele> getAlleles() { return alleles; }

    /**
     * @return the ploidy of this genotype
     */
    public int getPloidy() { return alleles.size(); }

    /**
     * @return true if all observed alleles are the same (regardless of whether they are ref or alt)
     */
    public boolean isHom() {

        // do not assume ploidy == 2

        if ( alleles.size() == 0 )
            return true;

        Allele first = alleles.get(0);
        for ( int i = 1; i < alleles.size(); i++ ) {
            if ( !first.equals(alleles.get(i)) )
                return false;
        }

        return true;
    }

    /**
     * @return true if we're het (observed alleles differ)
     */
    public boolean isHet() { return !isHom(); }

    /**
     * @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF)
     */
    public boolean isNoCall() {
        // TODO -- implement me
        return false;
    }

    /**
     * @return true if all alleles for this genotype are SNPs or reference
     */
//    public boolean isPointGenotype() {
////        for ( Allele allele : alleles ) {
////            if ( allele.isVariant() && !allele.isSNP() )
////                return false;
////        }
//        return true;
//    }

    /**
     * @return true if this is a variant genotype, false if it's reference
     */
    public boolean isVariant() {
        for ( Allele allele : alleles ) {
            if ( allele.isNonReference() )
                return true;
        }
        return false;
    }

    /**
     * @return the -1 * log10-based error estimate
     */
    public double getNegLog10PError() { return negLog10PError; }

    /**
     * @return the sample name
     */
    public String getSample() { return sample; }

    /**
     * Sets the sample name
     *
     * @param sample    the sample name
     */
    public void setSample(String sample) {
        this.sample = sample;
    }

    /**
     * @param ref the reference base
     *
     * @return this genotype as a variation
     */
    // TODO -- implement me
    // public Variation toVariation(char ref);
}
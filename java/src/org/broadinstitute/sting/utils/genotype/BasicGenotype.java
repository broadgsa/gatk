package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class BasicGenotype
 *         <p/>
 *         represents a basic genotype object.  That means that is
 *         an implementation for a basic genotype call, given the genotype
 *         string, the ref base, the confidence score, and the location. This
 *         class currently only represents point genotypes, not indels
 */
public class BasicGenotype implements Genotype {
    // the genotype string
    private String mGenotype;

    // our location
    private GenomeLoc mLocation;

    // the reference base.
    private char mRef;

    // the confidence score
    private double mNegLog10PError;

    /**
     * create a basic genotype, given the following fields
     *
     * @param location       the genomic location
     * @param genotype       the genotype, as a string, where ploidy = string.length
     * @param ref            the reference base as a char
     * @param negLog10PError the confidence score
     */
    public BasicGenotype(GenomeLoc location, String genotype, char ref, double negLog10PError) {
        mNegLog10PError = negLog10PError;
        mGenotype = genotype;
        mLocation = location;
        mRef = ref;
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the negitive log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return mNegLog10PError;
    }

    /**
     * get the bases that represent this genotype
     *
     * @return the bases, as a string
     */
    @Override
    public String getBases() {
        return mGenotype;
    }

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    @Override
    public int getPloidy() {
        return mGenotype.length();
    }

    /**
     * Returns true if both observed allele bases are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    @Override
    public boolean isHom() {
        if (mGenotype.length() < 1)
            return false;
        char base = mGenotype.charAt(0);
        for (char cur : mGenotype.toCharArray()) {
            if (base != cur) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns true if observed allele bases differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    @Override
    public boolean isHet() {
        if (mGenotype.length() < 1)
            return false;
        return !isHom();
    }

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    @Override
    public GenomeLoc getLocation() {
        return mLocation;
    }

    // set the location
    public void setLocation(GenomeLoc loc) {
        mLocation = loc;
    }

    /**
     * returns true if the genotype is a point genotype, false if it's a indel / deletion
     *
     * @return true is a SNP
     */
    @Override
    public boolean isPointGenotype() {
        return true;
    }

    /**
     * given the reference, are we a variant? (non-ref)
     *
     * @param ref the reference base or bases
     *
     * @return true if we're a variant
     */
    @Override
    public boolean isVariant(char ref) {
        return !(mGenotype.charAt(0) == ref && isHom());
    }

    /**
     * get the reference base.
     *
     * @return a character, representing the reference base
     */
    @Override
    public char getReference() {
        return mRef;
    }

    // set the reference base
    public void setReference(char refBase) {
        mRef = Character.toUpperCase(refBase);
    }

    /**
     * return this genotype as a variant
     *
     * @return the variant
     */
    @Override
    public Variation toVariation() {
        if (!isVariant(this.mRef)) throw new IllegalStateException("this genotype is not a variant");
        return new BasicVariation(this.getBases(), String.valueOf(mRef), this.getBases().length(), mLocation, mNegLog10PError);
    }

    /**
     * Turn a list of alleles into a genotype 
     * @param alleles the list of alleles
     * @return a string representation of this list
     */
    public static String alleleListToString(List<String> alleles) {
        StringBuilder builder = new StringBuilder();
        for (String allele : alleles)
            builder.append(allele);
        return builder.toString();
    }

}

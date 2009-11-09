package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Aug 4, 2009
 * Time: 6:46:09 PM
 * To change this template use File | Settings | File Templates.
 */
public enum DiploidGenotype {
    AA ('A', 'A'),
    AC ('A', 'C'),
    AG ('A', 'G'),
    AT ('A', 'T'),
    CC ('C', 'C'),
    CG ('C', 'G'),
    CT ('C', 'T'),
    GG ('G', 'G'),
    GT ('G', 'T'),
    TT ('T', 'T');

    public char base1, base2;
    private DiploidGenotype(char base1, char base2) {
        this.base1 = base1;
        this.base2 = base2;
    }

    public boolean isHomRef(char r) {
        return isHom() && r == this.toString().charAt(0);
    }

    public boolean isHomVar(char r) {
        return isHom() && r != this.toString().charAt(0);
    }

    public boolean isHetRef(char r) {
        return Utils.countOccurrences(r, this.toString()) == 1;
    }

    public boolean isHom() {
        return ! isHet();
    }

    public boolean isHet() {
        switch (this) {
            case AA:
            case CC:
            case GG:
            case TT: return false;
            default: return true;
        }
    }

    /**
     * create a diploid genotype, given a character to make into a hom genotype
     * @param hom the character to turn into a hom genotype, i.e. if it is A, then returned will be AA
     * @return the diploid genotype
     */
    public static DiploidGenotype createHomGenotype(char hom) {
        return DiploidGenotype.valueOf((String.valueOf(hom) + String.valueOf(hom)).toUpperCase());     
    }

    /**
     * get the genotype, given a string of 2 chars which may not necessarily be ordered correctly
     * @param genotype the string representation
     * @return the diploid genotype
     */
    public static DiploidGenotype unorderedValueOf(String genotype) {
        if ( genotype == null || genotype.length() != 2 )
            throw new IllegalArgumentException("Diploid genotypes are represented by 2 characters");
        if ( genotype.charAt(0) > genotype.charAt(1) )
            genotype = String.format("%c%c", genotype.charAt(1), genotype.charAt(0));
        return valueOf(genotype);
    }
}
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
        return isHom() && r == base1;
    }

    public boolean isHomVar(char r) {
        return isHom() && r != base1;
    }

    public boolean isHetRef(char r) {
        return Utils.countOccurrences(r, this.toString()) == 1;
    }

    public boolean isHom() {
        return ! isHet();
    }

    public boolean isHet() {
        return base1 != base2;
    }

    /**
     * create a diploid genotype, given a character to make into a hom genotype
     * @param hom the character to turn into a hom genotype, i.e. if it is A, then returned will be AA
     * @return the diploid genotype
     */
    public static DiploidGenotype createHomGenotype(char hom) {
        hom = Character.toUpperCase(hom);
        switch (hom) {
            case 'A': return DiploidGenotype.AA;
            case 'C': return DiploidGenotype.CC;
            case 'G': return DiploidGenotype.GG;
            case 'T': return DiploidGenotype.TT;
        }
        throw new IllegalArgumentException(hom + " is not a valid base character");
    }

    /**
     * get the genotype, given a string of 2 chars which may not necessarily be ordered correctly
     * @param base1 base1
     * @param base2 base2
     * @return the diploid genotype
     */
    public static DiploidGenotype unorderedValueOf(char base1, char base2) {
        if ( base1 > base2 ) {
            char temp = base1;
            base1 = base2;
            base2 = temp;
        }
        return valueOf(String.format("%c%c", base1, base2));
    }
}
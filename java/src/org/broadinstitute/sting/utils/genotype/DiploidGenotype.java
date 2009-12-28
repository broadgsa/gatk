package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;

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
        int index = BaseUtils.simpleBaseToBaseIndex(hom);
        if ( index == -1 )
            throw new IllegalArgumentException(hom + " is not a valid base character");
        return conversionMatrix[index][index];
    }

    /**
     * create a diploid genotype, given 2 chars which may not necessarily be ordered correctly
     * @param base1 base1
     * @param base2 base2
     * @return the diploid genotype
     */
    public static DiploidGenotype createDiploidGenotype(char base1, char base2) {
        int index1 = BaseUtils.simpleBaseToBaseIndex(base1);
        if ( index1 == -1 )
            throw new IllegalArgumentException(base1 + " is not a valid base character");
        int index2 = BaseUtils.simpleBaseToBaseIndex(base2);
        if ( index2 == -1 )
            throw new IllegalArgumentException(base2 + " is not a valid base character");
        return conversionMatrix[index1][index2];
    }

    private static final DiploidGenotype[][] conversionMatrix = {
            { DiploidGenotype.AA, DiploidGenotype.AC, DiploidGenotype.AG, DiploidGenotype.AT },
            { DiploidGenotype.AC, DiploidGenotype.CC, DiploidGenotype.CG, DiploidGenotype.CT },
            { DiploidGenotype.AG, DiploidGenotype.CG, DiploidGenotype.GG, DiploidGenotype.GT },
            { DiploidGenotype.AT, DiploidGenotype.CT, DiploidGenotype.GT, DiploidGenotype.TT }
    };
}
/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

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

    public byte base1, base2;

    @Deprecated
    private DiploidGenotype(char base1, char base2) {
        this((byte)base1, (byte)base2);
    }

    private DiploidGenotype(byte base1, byte base2) {
        this.base1 = base1;
        this.base2 = base2;
    }

    public boolean isHomRef(byte r) {
        return isHom() && r == base1;
    }

    public boolean isHomVar(byte r) {
        return isHom() && r != base1;
    }

    public boolean isHetRef(byte r) {
        if ( base1 == r )
            return r != base2;
        else
            return base2 == r;
        //return MathUtils.countOccurrences(r, this.toString()) == 1;
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
    public static DiploidGenotype createHomGenotype(byte hom) {
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
    public static DiploidGenotype createDiploidGenotype(byte base1, byte base2) {
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
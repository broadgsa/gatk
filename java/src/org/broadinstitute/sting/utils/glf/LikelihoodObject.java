package org.broadinstitute.sting.utils.glf;

import java.util.HashMap;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *
 * <p/>
 * Class LikelyhoodObject
 * <p/>
 * An object used to store likelyhood information for genotypes.  Genotype
 * likelihoods are assumed to be zero, unless set.  This allows the consumer
 * to make an empty LikelihoodObject, and just set those values which have
 * associated likelihood values.
 */
public class LikelihoodObject {

    // our possible genotypes, in order according to GLFv3
    public enum GENOTYPE {
        AA, AT, AC, AG, CC, CT, CG, GG, GT, TT
    }

    // how many genotypes we're storing
    public static final int genoTypeCount = GENOTYPE.values().length;

    // the associated negitive log likelihoods for each genotype
    final protected HashMap<GENOTYPE, Integer> likelihood = new HashMap<GENOTYPE, Integer>();

    /** create a blank likelihood object */
    public LikelihoodObject() {
        for (GENOTYPE type : GENOTYPE.values()) {
            likelihood.put(type, 255);
        }
    }

    // since the default case it to start with all infinity, we can choose any random base
    protected GENOTYPE greatest = GENOTYPE.AA;


    /**
     * create a likelyhood object, given an array of genotype scores in GLFv3 ordering
     * @param values an array of int's from 0 to 255, representing the negitive log likelihoods.
     */
    public LikelihoodObject(int[] values) {
        if (values.length != GENOTYPE.values().length) {
            throw new IllegalArgumentException("invalid array passed to LikelihoodObject, should be size " + GENOTYPE.values().length);
        }
        int index = 0;
        int lowestScore = 256;
        for (GENOTYPE type : GENOTYPE.values()) {
            if (values[index] < 0 || values[index] > 255) {
                throw new IllegalArgumentException("likelihood values must be greater or equal 0, and less then 256, value given: " + values[index]);
            }
            likelihood.put(type, values[index]);
            if (values[index] < lowestScore) {
                lowestScore = values[index];
                greatest = type;
            }
            ++index;
        }
    }

    /**
     * set the likelihood, given it's probability and the genotype
     * @param type the genotype
     * @param lh the likelihood as a double between 0 and 1, which is converted to a byte
     */
    public void setLikelihood(GENOTYPE type, int lh) {
        if (lh < 0 || lh > 255) {
            throw new IllegalArgumentException("supplied likelihood must be between 0 and 255");
        }
        likelihood.put(type,lh);
        if (lh < likelihood.get(this.greatest)) {
            this.greatest = type;
        }
    }

    /**
     * find the minimum likelihood value stored in the set.  This represents the most likely genotype,
     * since genotypes are represented as negitive log likeihoods
     * @return
     */
    public int getMinimumValue() {
        return likelihood.get(this.greatest);
    }

    /**
     * return a byte array representation of the likelihood object, in GLFv3 specified order.
     * The return type is short[] instead of byte[], since signed bytes only store -127 to 127,
     * not the 255 range we need.
     * @return a byte array of the genotype values
     */
    public short[] toByteArray() {
        short ret[] = new short[GENOTYPE.values().length];
        int index = 0;
        for (GENOTYPE type : GENOTYPE.values()) {
            ret[index] = (short)likelihood.get(type).intValue();
            ++index;
        }
        return ret;
    }


}

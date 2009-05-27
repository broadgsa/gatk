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

    // the associated probabilities for each genotype
    final protected HashMap<GENOTYPE, Integer> likelihood = new HashMap<GENOTYPE, Integer>();

    /** create a blank likelihood object */
    public LikelihoodObject() {
        for (GENOTYPE type : GENOTYPE.values()) {
            likelihood.put(type, 0);
        }
    }

    /**
     * create a likelyhood object, given an array of genotype scores in GLFv3 ordering
     * @param values
     */
    public LikelihoodObject(int[] values) {
        if (values.length != GENOTYPE.values().length) {
            throw new IllegalArgumentException("invalid array passed to LikelihoodObject, should be size " + GENOTYPE.values().length);
        }
        int index = 0;
        for (GENOTYPE type : GENOTYPE.values()) {
            if (values[index] < 0 || values[index] > 255) {
                throw new IllegalArgumentException("likelihood values must be greater or equal 0, and less then 256, value given: " + values[index]);
            }
            likelihood.put(type, values[index]);
            ++index;
        }
    }

    /**
     * set the likelihood, given it's probability and the genotype
     * @param type the genotype
     * @param likelyhood the likelihood as a float between 0 and 1, which is converted to a byte
     */
    public void setLikelihood(GENOTYPE type, float likelyhood) {
        likelihood.put(type,(int)Math.round(likelyhood*255.0));
    }

    /**
     * return a byte array representation of the likelihood object, in GLFv3 specified order
     * @return a byte array of the genotype values
     */
    public int[] toByteArray() {
        int ret[] = new int[GENOTYPE.values().length];
        int index = 0;
        for (GENOTYPE type : GENOTYPE.values()) {
            ret[index] = likelihood.get(type);
            ++index;
        }
        return ret;
    }


}

/*
 * Copyright (c) 2010.
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

import org.broadinstitute.sting.utils.variantcontext.Allele;

public class BiallelicGenotypeLikelihoods {

    private String sample;
    private double[] GLs;
    private Allele A, B;
    private int depth;

    /**
     * Create a new object for sample with given alleles and genotype likelihoods
     *
     * @param sample              sample name
     * @param A                   allele A
     * @param B                   allele B
     * @param log10AALikelihoods  AA likelihoods
     * @param log10ABLikelihoods  AB likelihoods
     * @param log10BBLikelihoods  BB likelihoods
     * @param depth               the read depth used in creating the likelihoods
     */
    public BiallelicGenotypeLikelihoods(String sample,
                                        Allele A,
                                        Allele B,
                                        double log10AALikelihoods,
                                        double log10ABLikelihoods,
                                        double log10BBLikelihoods,
                                        int depth) {
        this.sample = sample;
        this.A = A;
        this.B = B;
        this.GLs = new double[]{log10AALikelihoods, log10ABLikelihoods, log10BBLikelihoods};
        this.depth = depth;
    }

    public String getSample() {
        return sample;
    }

    public double getAALikelihoods() {
        return GLs[0];
    }

    public double getABLikelihoods() {
        return GLs[1];
    }

    public double getBBLikelihoods() {
        return GLs[2];
    }

    public double[] getLikelihoods() {
        return GLs;
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }

    public int getDepth() {
        return depth;
    }
}


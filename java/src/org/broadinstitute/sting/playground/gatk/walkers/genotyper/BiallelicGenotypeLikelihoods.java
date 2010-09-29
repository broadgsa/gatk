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

package org.broadinstitute.sting.playground.gatk.walkers.genotyper;

import org.broad.tribble.util.variantcontext.Allele;

public class BiallelicGenotypeLikelihoods {

    private String sample;
    private double[] GLs;
    private double[] GPs;
    private Allele A, B;

    /**
     * Create a new object for sample with given alleles and genotype likelihoods
     *
     * @param sample              sample name
     * @param A                   allele A
     * @param B                   allele B
     * @param log10AALikelihoods  AA likelihoods
     * @param log10ABLikelihoods  AB likelihoods
     * @param log10BBLikelihoods  BB likelihoods
     */
    public BiallelicGenotypeLikelihoods(String sample,
                                        Allele A,
                                        Allele B,
                                        double log10AALikelihoods,
                                        double log10ABLikelihoods,
                                        double log10BBLikelihoods) {
        this.sample = sample;
        this.A = A;
        this.B = B;
        this.GLs = new double[]{log10AALikelihoods, log10ABLikelihoods, log10BBLikelihoods};
        this.GPs = new double[]{log10AALikelihoods, log10ABLikelihoods, log10BBLikelihoods};
    }

    /**
     * Create a new object for sample with given alleles and genotype likelihoods & posteriors
     *
     * @param sample              sample name
     * @param A                   allele A
     * @param B                   allele B
     * @param log10AALikelihoods  AA likelihoods
     * @param log10ABLikelihoods  AB likelihoods
     * @param log10BBLikelihoods  BB likelihoods
     * @param log10AAPosteriors   AA posteriors
     * @param log10ABPosteriors   AB posteriors
     * @param log10BBPosteriors   BB posteriors
     */
    public BiallelicGenotypeLikelihoods(String sample,
                                        Allele A,
                                        Allele B,
                                        double log10AALikelihoods,
                                        double log10ABLikelihoods,
                                        double log10BBLikelihoods,
                                        double log10AAPosteriors,
                                        double log10ABPosteriors,
                                        double log10BBPosteriors) {
        this.sample = sample;
        this.A = A;
        this.B = B;
        this.GLs = new double[]{log10AALikelihoods, log10ABLikelihoods, log10BBLikelihoods};
        this.GPs = new double[]{log10AAPosteriors, log10ABPosteriors, log10BBPosteriors};
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

    public double getAAPosteriors() {
        return GPs[0];
    }

    public double getABPosteriors() {
        return GPs[1];
    }

    public double getBBPosteriors() {
        return GPs[2];
    }

    public double[] getPosteriors() {
        return GPs;
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }
}

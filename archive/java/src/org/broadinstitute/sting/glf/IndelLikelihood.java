package org.broadinstitute.sting.utils.genotype;


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
 *         <p/>
 *         Class IndelLikelihood
 *         <p/>
 *         The representation of an indel allele.
 */
// TODO -- DELETE ME GLF
public class IndelLikelihood {

    protected double loglikelihood;
    protected short lengthOfIndel;
    protected short[] indelSequence;

    /**
     * Create a likelihood object for an indel call
     *
     * @param likelihood    the likelihood (represented as a negitive log likelihood,
     *                      with a ceiling of 255.
     * @param indelSequence the indel sequence, not null terminated
     */
    public IndelLikelihood( byte likelihood, String indelSequence ) {
        this.loglikelihood = likelihood;
        this.lengthOfIndel = (short)indelSequence.length();
        this.indelSequence = new short[indelSequence.length()];
        for (int tmp = 0; tmp < indelSequence.length(); tmp++) {
            this.indelSequence[tmp] = (short)indelSequence.charAt(tmp);    
        }
    }

    /**
     * getter methods
     */
    
    public double getLikelihood() {
        return loglikelihood;
    }

    public short getLengthOfIndel() {
        return lengthOfIndel;
    }

    public short[] getIndelSequence() {
        return indelSequence;
    }
}

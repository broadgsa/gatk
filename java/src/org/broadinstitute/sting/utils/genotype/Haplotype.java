/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.genotype;

import java.util.Arrays;

public class Haplotype {
    protected byte[] bases = null;
    protected double[] quals = null;

    /**
     * Create a simple consensus sequence with provided bases and a uniform quality over all bases of qual
     *
     * @param bases bases
     * @param qual  qual
     */
    public Haplotype(byte[] bases, int qual) {
        this.bases = bases;
        quals = new double[bases.length];
        Arrays.fill(quals, (double)qual);
    }

    public Haplotype(byte[] bases, double[] quals) {
        this.bases = bases;
        this.quals = quals;
    }

    public Haplotype(String bases, double[] quals) {
        this.bases = bases.getBytes();
        this.quals = quals;
    }

    public Haplotype(byte[] bases) {
        this(bases, 0);
    }


    public String toString() { return new String(this.bases); }

    public double getQualitySum() {
        double s = 0;
        for (int k=0; k < bases.length; k++) {
            s += quals[k];
        }
        return s;
    }

    public double[] getQuals() {
        return quals;
    }
    public byte[] getBasesAsBytes() {
        return bases;
    }
}

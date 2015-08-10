/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.pairhmm;

import com.google.java.contract.Requires;

/**
 * Superclass for PairHMM that want to use a full read x haplotype matrix for their match, insertion, and deletion matrix
 *
 * User: rpoplin
 * Date: 10/16/12
 */
abstract class N2MemoryPairHMM extends PairHMM {
    protected double[][] transition = null; // The transition probabilities cache
    protected double[][] prior = null;      // The prior probabilities cache
    protected double[][] matchMatrix = null;
    protected double[][] insertionMatrix = null;
    protected double[][] deletionMatrix = null;

    // only used for debugging purposes
    protected boolean doNotUseTristateCorrection = false;

    public void doNotUseTristateCorrection() {
        doNotUseTristateCorrection = true;
    }

    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     *
     * Note: Do not worry about padding, just provide the true max length of the read and haplotype. The HMM will take care of the padding.
     *
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     */
    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        matchMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
        insertionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
        deletionMatrix = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];

        transition = PairHMMModel.createTransitionMatrix(maxReadLength);
        prior = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];
    }

    /**
     * Print out the core hmm matrices for debugging
     */
    protected void dumpMatrices() {
        dumpMatrix("matchMetricArray", matchMatrix);
        dumpMatrix("insertionMatrix", insertionMatrix);
        dumpMatrix("deletionMatrix", deletionMatrix);
    }

    /**
     * Print out in a human readable form the matrix for debugging
     * @param name the name of this matrix
     * @param matrix the matrix of values
     */
    @Requires({"name != null", "matrix != null"})
    private void dumpMatrix(final String name, final double[][] matrix) {
        System.out.printf("%s%n", name);
        for ( int i = 0; i < matrix.length; i++) {
            System.out.printf("\t%s[%d]", name, i);
            for ( int j = 0; j < matrix[i].length; j++ ) {
                if ( Double.isInfinite(matrix[i][j]) )
                    System.out.printf(" %15s", String.format("%f", matrix[i][j]));
                else
                    System.out.printf(" % 15.5e", matrix[i][j]);
            }
            System.out.println();
        }
    }
}

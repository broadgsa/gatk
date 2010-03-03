package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 24, 2010
 */

public class VariantTreeNode {

    public VariantTreeNode left;
    public VariantTreeNode right;
    public int cutDim;
    public double cutValue;
    public VariantDatum[] variants;

    private final int minBinSize = 8000; // BUGBUG: must be larger than number of kNN

    public VariantTreeNode() {
        left = null;
        right = null;
        variants = null;
        cutDim = -1;
        cutValue = -1;
    }

    public void cutData( final VariantDatum[] data, final int depth, final int lastCutDepth, final int numAnnotations ) {

        cutDim = depth % numAnnotations;

        if( depth != lastCutDepth && (cutDim == (lastCutDepth % numAnnotations)) ) { // Base case: we've tried to cut on all the annotations
            variants = data;
            return;
        }

        final double[] values = new double[data.length];
        for( int iii = 0; iii < data.length; iii++ ) {
            values[iii] = data[iii].annotations[cutDim];
        }

        final double[] sortedValues = values.clone();
        VariantTree.quickSort( sortedValues, 0, values.length-1 ); // BUGBUG: values.length or values.length-1

        final int lowPivotIndex = Math.round(0.40f * sortedValues.length);
        final int highPivotIndex = Math.round(0.60f * sortedValues.length);
        final double lowPivot = sortedValues[lowPivotIndex];
        final double highPivot = sortedValues[highPivotIndex];
        cutValue = highPivot;

        int numLow = 0;
        int numHigh = 0;
        for( int iii = 0; iii < data.length; iii++ ) {
            if( values[iii] < highPivot ) { numLow++; }
            if( values[iii] >= lowPivot ) { numHigh++; }
        }

        // If cutting here makes the bin too small then don't cut
        if( numLow < minBinSize || numHigh < minBinSize  || (numLow == numHigh && numLow == data.length) ) {
            cutValue = sortedValues[0];
            right = new VariantTreeNode();
            right.cutData(data, depth+1, lastCutDepth, numAnnotations);
        } else {
            final VariantDatum[] leftData = new VariantDatum[numLow];
            final VariantDatum[] rightData = new VariantDatum[numHigh];
            int leftIndex = 0;
            int rightIndex = 0;
            for( int iii = 0; iii < data.length; iii++ ) {
                if( values[iii] < highPivot ) { leftData[leftIndex++] = data[iii]; }
                if( values[iii] >= lowPivot ) { rightData[rightIndex++] = data[iii]; }
            }

            left = new VariantTreeNode();
            right = new VariantTreeNode();
            left.cutData(leftData, depth+1, depth, numAnnotations);
            right.cutData(rightData, depth+1, depth, numAnnotations);
        }
    }

}

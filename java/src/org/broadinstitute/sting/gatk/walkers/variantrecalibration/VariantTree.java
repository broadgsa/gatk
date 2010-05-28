package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.utils.StingException;

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

public class VariantTree {
    private final VariantTreeNode root;
    private final int numKNN;

    public VariantTree( final int _numKNN ) {
        root = new VariantTreeNode();
        numKNN = _numKNN;
    }

    public void createTreeFromData( VariantDatum[] data ) {
        root.cutData( data, 0, 0, data[0].annotations.length );
    }

    public double calcNeighborhoodTITV( final VariantDatum variant ) {

        double[] distances;

        // Grab the subset of points that are approximately near this point
        final VariantDatum[] data = getBin( variant.annotations, root );
        if( data.length < numKNN ) {
            throw new StingException( "Bin is too small. Should be > " + numKNN );
        }

        // Find the X nearest points in the subset
        final double[] originalDistances = calcDistances( variant.annotations, data );
        distances = originalDistances.clone();
        quickSort( distances, 0, distances.length-1 ); // BUGBUG: distances.length or distances.length-1

        final double minDistance = distances[numKNN - 1];

        // Calculate probability of being true based on this set of SNPs
        int numTi = 0;
        int numTv = 0;
        for( int iii = 0; iii < distances.length; iii++ ) {
            if( originalDistances[iii] <= minDistance ) {
                if( data[iii].isTransition ) { numTi++; }
                else { numTv++; }
            }
        }

        return ((double) numTi) / ((double) numTv);
    }

    private VariantDatum[] getBin( final double[] variant, final VariantTreeNode node ) {
        if( node.variants != null ) {
            return node.variants;
        } else {
            if( variant[node.cutDim] < node.cutValue ) {
                return getBin( variant, node.left );
            } else {
                return getBin( variant, node.right );
            }
        }
    }

    private double[] calcDistances( final double[] variant, final VariantDatum[] data ) {
        final double[] distSquared = new double[data.length];
        int iii = 0;
        for( final VariantDatum variantDatum : data ) {
            distSquared[iii] = 0.0;
            int jjj = 0;
            for( final double value : variantDatum.annotations) {
                final double diff = variant[jjj] - value;
                distSquared[iii] += ( diff * diff );
                jjj++;
            }
            iii++;
        }

        return distSquared;
    }

    public static int partition(final double arr[], final int left, final int right)
    {
          int i = left, j = right;
          double tmp;
          final double pivot = arr[(left + right) / 2];

          while (i <= j) {
                while (arr[i] < pivot) {
                      i++;
                }
                while (arr[j] > pivot) {
                      j--;
                }
                if (i <= j) {
                      tmp = arr[i];
                      arr[i] = arr[j];
                      arr[j] = tmp;
                      i++;
                      j--;
                }
          }

          return i;
    }

    public static void quickSort(final double arr[], final int left, final int right) {
          final int index = partition(arr, left, right);
          if (left < index - 1) {
                quickSort(arr, left, index - 1);
          }
          if (index < right) {
                quickSort(arr, index, right);
          }
    }


}

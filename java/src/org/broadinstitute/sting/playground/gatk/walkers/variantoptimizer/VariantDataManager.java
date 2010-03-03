package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.utils.ExpandingArrayList;
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
 * Date: Feb 26, 2010
 */

public class VariantDataManager {
    public final VariantDatum[] data;
    public final int numVariants;
    public final int numAnnotations;
    public final double[] varianceVector;
    public boolean isNormalized;
    private final ExpandingArrayList<String> annotationKeys;

    public VariantDataManager( final ExpandingArrayList<VariantDatum> dataList, final ExpandingArrayList<String> _annotationKeys ) {
        numVariants = dataList.size();
        data = dataList.toArray( new VariantDatum[numVariants]);
        if( numVariants <= 0 ) {
            throw new StingException( "There are zero variants! (or possibly a problem with integer overflow)" );
        }
        numAnnotations = data[0].annotations.length;
        if( numAnnotations <= 0 ) {
            throw new StingException( "There are zero annotations! (or possibly a problem with integer overflow)" );
        }
        varianceVector = new double[numAnnotations];
        isNormalized = false;
        annotationKeys = _annotationKeys;
    }

    public void normalizeData() {
        for( int iii = 0; iii < numAnnotations; iii++ ) {
            final double theMean = mean(data, iii);
            final double theSTD = standardDeviation(data, theMean, iii);
            System.out.println( (iii == numAnnotations-1 ? "QUAL" : annotationKeys.get(iii)) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );
            varianceVector[iii] = theSTD * theSTD;
            for( int jjj=0; jjj<numVariants; jjj++ ) {
                data[jjj].annotations[iii] = ( data[jjj].annotations[iii] - theMean ) / theSTD;
            }
        }
        isNormalized = true; // Each data point is now [ (x - mean) / standard deviation ]
    }

    private static double mean( final VariantDatum[] data, final int index ) {
        double sum = 0.0;
        final int numVars = data.length;
        for( int iii = 0; iii < numVars; iii++ ) {
            sum += (data[iii].annotations[index] / ((double) numVars));
        }
        return sum;
    }

    private static double standardDeviation( final VariantDatum[] data, final double mean, final int index ) {
        double sum = 0.0;
        final int numVars = data.length;
        for( int iii = 0; iii < numVars; iii++ ) {
            sum += ( ((data[iii].annotations[index] - mean)*(data[iii].annotations[index] - mean)) / ((double) numVars));
        }
        return Math.sqrt( sum );
    }

}

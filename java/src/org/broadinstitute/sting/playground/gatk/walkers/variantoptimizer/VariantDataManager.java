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

    public VariantDataManager( final ExpandingArrayList<VariantDatum> dataList ) {
        numVariants = dataList.size();
        data = dataList.toArray( new VariantDatum[numVariants]);
        if( numVariants <= 0 ) {
            throw new StingException( "There are zero variants! (or possibly problem with integer overflow)" );
        }
        numAnnotations = data[0].annotations.length;
        if( numAnnotations <= 0 ) {
            throw new StingException( "There are zero annotations! (or possibly problem with integer overflow)" );
        }
        varianceVector = new double[numAnnotations];
        isNormalized = false;
    }

    public void normalizeData() {
        for( int iii = 0; iii < numAnnotations; iii++ ) {
            final double theMean = mean(data, iii);
            final double theVariance = variance(data, theMean, iii);
            varianceVector[iii] = theVariance * theVariance; // BUGBUG: model's use this as sigma^2 instead of sigma
            for( int jjj=0; jjj<numVariants; jjj++ ) {
                data[jjj].annotations[iii] = ( data[jjj].annotations[iii] - theMean ) / theVariance;
            }
        }
        isNormalized = true;
    }

    private static double mean( final VariantDatum[] data, final int index ) {
        double sum = 0.0f;
        final int numVars = data.length;
        for( int iii = 0; iii < numVars; iii++ ) {
            sum += (data[iii].annotations[index] / ((double) numVars));
        }
        return sum;
    }

    private static double variance( final VariantDatum[] data, final double mean, final int index ) {
        double sum = 0.0f;
        final int numVars = data.length;
        for( int iii = 0; iii < numVars; iii++ ) {
            sum += ( ((data[iii].annotations[index] - mean)*(data[iii].annotations[index] - mean)) / ((double) numVars));
        }
        return Math.sqrt( sum );
    }

}

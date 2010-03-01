package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.StingException;

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

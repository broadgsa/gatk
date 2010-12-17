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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.Utils;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 26, 2010
 */

public class VariantDataManager {
    public enum TailType {
        TWO_TAILED,
        SMALL_IS_GOOD,
        BIG_IS_GOOD
    }

    protected final static Logger logger = Logger.getLogger(VariantDataManager.class);

    public final VariantDatum[] data;
    public final int numVariants;
    public final int numAnnotations;
    public final double[] meanVector;
    public final double[] varianceVector; // This is really the standard deviation
    public final TailType[] tailTypes;
    public boolean isNormalized;
    public final ExpandingArrayList<String> annotationKeys;

    public VariantDataManager( final ExpandingArrayList<VariantDatum> dataList, final ExpandingArrayList<String> _annotationKeys ) {
        numVariants = dataList.size();
        data = dataList.toArray( new VariantDatum[numVariants] );
        if( numVariants <= 0 ) {
            throw new UserException.BadInput("There are zero variants with > 0 clustering weight! Please provide sets of known polymorphic loci to be used as training data using the dbsnp, hapmap, or 1kg rod bindings. Clustering weights can be specified using -weightDBSNP, -weightHapMap, and -weight1KG" );
        }
        if( _annotationKeys == null ) {
            numAnnotations = 0;
            meanVector = null;
            varianceVector = null;
            tailTypes = null;
        } else {
            numAnnotations = _annotationKeys.size();
            if( numAnnotations <= 0 ) {
                throw new UserException.BadInput("There are zero annotations.  At least one annotation must be provided to use this walker!" );
            }
            meanVector = new double[numAnnotations];
            varianceVector = new double[numAnnotations];
            tailTypes = new TailType[numAnnotations];
            Arrays.fill(tailTypes, TailType.TWO_TAILED);
            isNormalized = false;
        }
        annotationKeys = _annotationKeys;
    }

    public VariantDataManager( final ExpandingArrayList<String> annotationLines ) {
        data = null;
        numVariants = 0;
        numAnnotations = annotationLines.size();
        meanVector = new double[numAnnotations];
        varianceVector = new double[numAnnotations];
        isNormalized = true;
        annotationKeys = new ExpandingArrayList<String>();
        tailTypes = new TailType[numAnnotations];
        Arrays.fill(tailTypes, TailType.TWO_TAILED);

        int jjj = 0;
        for( final String line : annotationLines ) {
            final String[] vals = line.split(",");
            annotationKeys.add(vals[1]);
            meanVector[jjj] = Double.parseDouble(vals[2]);
            varianceVector[jjj] = Double.parseDouble(vals[3]);
            jjj++;
        }
    }

    //
    // code for dealing with TailTypes
    //
    
    public void setAnnotationType(String annotation, TailType t ) {
        int i = annotationKeys.indexOf(annotation);
        if ( i == -1 ) throw new UserException.BadInput("Attempt to modify unknown annotation " + annotation + ": cluster file has " + Utils.join(",", annotationKeys));
        tailTypes[i] = t;
        logger.info("Updating annotation " + annotation + " to have tail type " + t);
    }

    public TailType getAnnotationType(int i) {
        return tailTypes[i];
    }


    public void normalizeData() {
        boolean foundZeroVarianceAnnotation = false;
        for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
            final double theMean = mean(data, jjj);
            final double theSTD = standardDeviation(data, theMean, jjj);
            logger.info( annotationKeys.get(jjj) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );
            if( theSTD < 1E-8 ) {
                foundZeroVarianceAnnotation = true;
                logger.warn("Zero variance is a problem: standard deviation = " + theSTD + " User must -exclude annotations with zero variance. Annotation = " + (jjj == numAnnotations-1 ? "QUAL" : annotationKeys.get(jjj)));
            } else if( theSTD < 1E-2 ) {
                logger.warn("Warning! Tiny variance. It is strongly recommended that you -exclude " + annotationKeys.get(jjj));
            }
            meanVector[jjj] = theMean;
            varianceVector[jjj] = theSTD;
            for( int iii = 0; iii < numVariants; iii++ ) {
                data[iii].annotations[jjj] = ( data[iii].annotations[jjj] - theMean ) / theSTD;
            }
        }
        isNormalized = true; // Each data point is now [ (x - mean) / standard deviation ]
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput("Found annotations with zero variance. They must be excluded before proceeding.");
        }
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

    public void printClusterFileHeader( PrintStream outputFile ) {
        for( int jjj = 0; jjj < numAnnotations; jjj++ ) {
            outputFile.println(String.format("@!ANNOTATION," + annotationKeys.get(jjj) + ",%.8f,%.8f", meanVector[jjj], varianceVector[jjj]));
        }
    }
}

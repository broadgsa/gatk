package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDataManager {
    public ExpandingArrayList<VariantDatum> data;
    public final double[] meanVector;
    public final double[] varianceVector; // This is really the standard deviation
    public final ArrayList<String> annotationKeys;

    private final static long RANDOM_SEED = 83409701;
    private final static Random rand = new Random( RANDOM_SEED ); // this is going to cause problems if it is ever used in an integration test

    protected final static Logger logger = Logger.getLogger(VariantDataManager.class);

    public VariantDataManager( final ExpandingArrayList<VariantDatum> data, final List<String> annotationKeys ) {
        this.data = data;
        this.annotationKeys = new ArrayList<String>( annotationKeys );
        meanVector = new double[this.annotationKeys.size()];
        varianceVector = new double[this.annotationKeys.size()];
    }

    public VariantDataManager( final List<String> annotationLines ) {
        data = null;
        annotationKeys = new ArrayList<String>( annotationLines.size() );
        meanVector = new double[annotationLines.size()];
        varianceVector = new double[annotationLines.size()];

        int jjj = 0;
        for( final String line : annotationLines ) {
            final String[] vals = line.split(",");
            annotationKeys.add(vals[1]);
            meanVector[jjj] = Double.parseDouble(vals[2]);
            varianceVector[jjj] = Double.parseDouble(vals[3]);
            jjj++;
        }
    }

    public void add( final VariantDatum datum ) {
        data.add( datum );
    }

    public void normalizeData() {
        boolean foundZeroVarianceAnnotation = false;
        for( int jjj = 0; jjj < meanVector.length; jjj++ ) { //BUGBUG: to clean up
            final double theMean = mean(jjj); //BUGBUG: to clean up
            final double theSTD = standardDeviation(theMean, jjj); //BUGBUG: to clean up
            logger.info( annotationKeys.get(jjj) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );
            if( theSTD < 1E-8 ) { //BUGBUG: output recreated here to match the integration tests
                foundZeroVarianceAnnotation = true;
            } else if( theSTD < 1E-2 ) {
                logger.warn("Warning! Tiny variance. It is strongly recommended that you -exclude " + annotationKeys.get(jjj));
            }
            meanVector[jjj] = theMean;
            varianceVector[jjj] = theSTD;
            for(final VariantDatum datum : data ) {
                datum.annotations[jjj] = ( datum.annotations[jjj] - theMean ) / theSTD; // Each data point is now [ (x - mean) / standard deviation ]
            }
        }
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput("Found annotations with zero variance. They must be excluded before proceeding.");
        }
    }

    private double mean( final int index ) {
        double sum = 0.0;
        final int numVars = data.size();
        for( final VariantDatum datum : data ) {
            sum += (datum.annotations[index] / ((double) numVars));
        }
        return sum;
    }

    private double standardDeviation( final double mean, final int index ) {
        double sum = 0.0;
        final int numVars = data.size();
        for( final VariantDatum datum : data ) {
            sum += ( ((datum.annotations[index] - mean)*(datum.annotations[index] - mean)) / ((double) numVars));
        }
        return Math.sqrt( sum );
    }

    public void trimDataBySTD( final double STD_THRESHOLD ) {
        final ExpandingArrayList<VariantDatum> dataToRemove = new ExpandingArrayList<VariantDatum>();
        for( final VariantDatum datum : data ) {
            boolean remove = false;
            for( double val : datum.annotations ) {
                if( Math.abs(val) > STD_THRESHOLD ) {
                    remove = true;
                }
            }
            if( remove ) { dataToRemove.add( datum ); }
        }
        data.removeAll( dataToRemove );
    }

    public void printClusterFileHeader( PrintStream outputFile ) {
        for( int jjj = 0; jjj < annotationKeys.size(); jjj++ ) {
            outputFile.println(String.format("@!ANNOTATION," + annotationKeys.get(jjj) + ",%.8f,%.8f", meanVector[jjj], varianceVector[jjj]));
        }
    }

    public static double decodeAnnotation( final GenomeLocParser genomeLocParser, final String annotationKey, final VariantContext vc, final boolean jitter ) {
        double value;
        if( jitter && annotationKey.equalsIgnoreCase("HRUN") ) { // HRun values must be jittered a bit to work in this GMM
            value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            value += -0.25 + 0.5 * rand.nextDouble(); //BUGBUG: recreated here to match the integration tests
        } else if( annotationKey.equals("QUAL") ) {
            value = vc.getPhredScaledQual();
        } else {
            try {
                value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            } catch( Exception e ) {
                throw new UserException.MalformedFile(vc.getSource(), "No double value detected for annotation = " + annotationKey + " in variant at " + VariantContextUtils.getLocation(genomeLocParser,vc) + ", reported annotation value = " + vc.getAttribute( annotationKey ), e );
            }
        }
        return value;
    }

    public double[] normalizeDatum( final double[] annotationValues ) {
        for(int iii = 0; iii < annotationValues.length; iii++) {
            annotationValues[iii] = (annotationValues[iii] - meanVector[iii]) / varianceVector[iii];
        }
        return annotationValues;
    }
}

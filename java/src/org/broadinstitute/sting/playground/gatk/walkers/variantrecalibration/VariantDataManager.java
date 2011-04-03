package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDataManager {
    private ExpandingArrayList<VariantDatum> data;
    private final double[] meanVector;
    private final double[] varianceVector; // this is really the standard deviation
    public final ArrayList<String> annotationKeys;
    private final ExpandingArrayList<TrainingSet> trainingSets;
    private final VariantRecalibratorArgumentCollection VRAC;

    private final static long RANDOM_SEED = 83409701;
    private final static Random rand = new Random( RANDOM_SEED ); // this is going to cause problems if it is ever used in an integration test, planning to get rid of HRun anyway

    protected final static Logger logger = Logger.getLogger(VariantDataManager.class);

    public VariantDataManager( final List<String> annotationKeys, final VariantRecalibratorArgumentCollection VRAC ) {
        this.data = null;
        this.annotationKeys = new ArrayList<String>( annotationKeys );
        this.VRAC = VRAC;
        meanVector = new double[this.annotationKeys.size()];
        varianceVector = new double[this.annotationKeys.size()];
        trainingSets = new ExpandingArrayList<TrainingSet>();
    }

    public void setData( final ExpandingArrayList<VariantDatum> data ) {
        this.data = data;
    }

    public ExpandingArrayList<VariantDatum> getData() {
        return data;
    }

    public void normalizeData() {
        boolean foundZeroVarianceAnnotation = false;
        for( int jjj = 0; jjj < meanVector.length; jjj++ ) { //BUGBUG: to clean up
            final double theMean = mean(jjj); //BUGBUG: to clean up
            final double theSTD = standardDeviation(theMean, jjj); //BUGBUG: to clean up
            logger.info( annotationKeys.get(jjj) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );
            foundZeroVarianceAnnotation = foundZeroVarianceAnnotation || (theSTD < 1E-8);
            meanVector[jjj] = theMean;
            varianceVector[jjj] = theSTD;
            for( final VariantDatum datum : data ) {
                datum.annotations[jjj] = ( datum.annotations[jjj] - theMean ) / theSTD; // Each data point is now [ (x - mean) / standard deviation ]
            }
        }
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput( "Found annotations with zero variance. They must be excluded before proceeding." );
        }
    }

    public void addTrainingSet( final TrainingSet trainingSet ) {
        trainingSets.add( trainingSet );
    }

    public boolean checkHasTrainingSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isTraining ) { return true; }
        }
        return false;
    }

    public boolean checkHasTruthSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isTruth ) { return true; }
        }
        return false;
    }

    public boolean checkHasKnownSet() {
        for( final TrainingSet trainingSet : trainingSets ) {
            if( trainingSet.isKnown ) { return true; }
        }
        return false;
    }

    public ExpandingArrayList<VariantDatum> getTrainingData() {
        final ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();
        for( final VariantDatum datum : data ) {
            if( datum.atTrainingSite && datum.originalQual > VRAC.QUAL_THRESHOLD ) {
                trainingData.add( datum );
            }
        }
        trimDataBySTD( trainingData, VRAC.STD_THRESHOLD );
        logger.info( "Training with " + trainingData.size() + " variants found in the training set(s)." );
        return trainingData;
    }

    public ExpandingArrayList<VariantDatum> selectWorstVariants( final double bottomPercentage ) {
        Collections.sort( data );
        final ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();
        trainingData.addAll( data.subList(0, Math.round((float)bottomPercentage * data.size())) );
        logger.info( "Training with worst " + bottomPercentage * 100.0f + "% of data --> " + trainingData.size() + " variants with LOD <= " + String.format("%.4f", data.get(Math.round((float)bottomPercentage * data.size())).lod) + "." );
        return trainingData;
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

    public static void trimDataBySTD( final ExpandingArrayList<VariantDatum> listData, final double STD_THRESHOLD ) {
        final ExpandingArrayList<VariantDatum> dataToRemove = new ExpandingArrayList<VariantDatum>();
        for( final VariantDatum datum : listData ) {
            boolean remove = false;
            for( final double val : datum.annotations ) {
                remove = remove || (Math.abs(val) > STD_THRESHOLD);
            }
            if( remove ) { dataToRemove.add( datum ); }
        }
        listData.removeAll( dataToRemove );
    }

    public double[] decodeAnnotations( final GenomeLocParser genomeLocParser, final VariantContext vc, final boolean jitter ) {
        final double[] annotations = new double[annotationKeys.size()];
        int iii = 0;
        for( final String key : annotationKeys ) {
            annotations[iii++] = decodeAnnotation( genomeLocParser, key, vc, jitter );
        }
        return annotations;
    }

    private static double decodeAnnotation( final GenomeLocParser genomeLocParser, final String annotationKey, final VariantContext vc, final boolean jitter ) {
        double value;
        if( jitter && annotationKey.equalsIgnoreCase("HRUN") ) { // HRun values must be jittered a bit to work in this GMM
            value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            value += -0.25 + 0.5 * rand.nextDouble();
        } else if( annotationKey.equals("QUAL") ) {
            value = vc.getPhredScaledQual();
        } else {
            try {
                value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            } catch( final Exception e ) {
                throw new UserException.MalformedFile( vc.getSource(), "No double value detected for annotation = " + annotationKey + " in variant at " + VariantContextUtils.getLocation(genomeLocParser,vc) + ", reported annotation value = " + vc.getAttribute( annotationKey ), e );
            }
        }
        return value;
    }

    public void parseTrainingSets( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context, final VariantContext evalVC, final VariantDatum datum, final boolean TRUST_ALL_POLYMORPHIC ) {
        datum.isKnown = false;
        datum.atTruthSite = false;
        datum.atTrainingSite = false;
        datum.prior = 3.0;
        for( final TrainingSet trainingSet : trainingSets ) {
            final Collection<VariantContext> vcs = tracker.getVariantContexts( ref, trainingSet.name, null, context.getLocation(), false, true );
            final VariantContext trainVC = ( vcs.size() != 0 ? vcs.iterator().next() : null );
            if( trainVC != null && trainVC.isVariant() && trainVC.isNotFiltered() && ((evalVC.isSNP() && trainVC.isSNP()) || (evalVC.isIndel() && trainVC.isIndel())) && (TRUST_ALL_POLYMORPHIC || !trainVC.hasGenotypes() || trainVC.isPolymorphic()) ) {
                datum.isKnown = datum.isKnown || trainingSet.isKnown;
                datum.atTruthSite = datum.atTruthSite || trainingSet.isTruth;
                datum.atTrainingSite = datum.atTrainingSite || trainingSet.isTraining;
                datum.prior = Math.max( datum.prior, trainingSet.prior );
            }
        }
    }

    public void writeOutRecalibrationTable( final PrintStream RECAL_FILE ) {
        for( final VariantDatum datum : data ) {
            RECAL_FILE.println(String.format("%s,%d,%d,%.4f", datum.pos.getContig(), datum.pos.getStart(), datum.pos.getStop(), datum.lod));
        }
    }
}

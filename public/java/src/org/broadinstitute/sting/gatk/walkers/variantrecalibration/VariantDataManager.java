/*
 * Copyright (c) 2011 The Broad Institute
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
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDataManager {
    private ExpandingArrayList<VariantDatum> data;
    private final double[] meanVector;
    private final double[] varianceVector; // this is really the standard deviation
    public final List<String> annotationKeys;
    private final VariantRecalibratorArgumentCollection VRAC;
    protected final static Logger logger = Logger.getLogger(VariantDataManager.class);
    protected final List<TrainingSet> trainingSets;

    public VariantDataManager( final List<String> annotationKeys, final VariantRecalibratorArgumentCollection VRAC ) {
        this.data = null;
        this.annotationKeys = new ArrayList<String>( annotationKeys );
        this.VRAC = VRAC;
        meanVector = new double[this.annotationKeys.size()];
        varianceVector = new double[this.annotationKeys.size()];
        trainingSets = new ArrayList<TrainingSet>();
    }

    public void setData( final ExpandingArrayList<VariantDatum> data ) {
        this.data = data;
    }

    public ExpandingArrayList<VariantDatum> getData() {
        return data;
    }

    public void normalizeData() {
        boolean foundZeroVarianceAnnotation = false;
        for( int iii = 0; iii < meanVector.length; iii++ ) {
            final double theMean = mean(iii);
            final double theSTD = standardDeviation(theMean, iii);
            logger.info( annotationKeys.get(iii) + String.format(": \t mean = %.2f\t standard deviation = %.2f", theMean, theSTD) );
            if( Double.isNaN(theMean) ) {
                throw new UserException.BadInput("Values for " + annotationKeys.get(iii) + " annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations. See http://www.broadinstitute.org/gsa/wiki/index.php/VariantAnnotator");
            }

            foundZeroVarianceAnnotation = foundZeroVarianceAnnotation || (theSTD < 1E-6);
            meanVector[iii] = theMean;
            varianceVector[iii] = theSTD;
            for( final VariantDatum datum : data ) {
                // Transform each data point via: (x - mean) / standard deviation
                datum.annotations[iii] = ( datum.isNull[iii] ? GenomeAnalysisEngine.getRandomGenerator().nextGaussian() : ( datum.annotations[iii] - theMean ) / theSTD );
            }
        }
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput( "Found annotations with zero variance. They must be excluded before proceeding." );
        }

        // trim data by standard deviation threshold and mark failing data for exclusion later
        for( final VariantDatum datum : data ) {
            boolean remove = false;
            for( final double val : datum.annotations ) {
                remove = remove || (Math.abs(val) > VRAC.STD_THRESHOLD);
            }
            datum.failingSTDThreshold = remove;
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
            if( datum.atTrainingSite && !datum.failingSTDThreshold && datum.originalQual > VRAC.QUAL_THRESHOLD ) {
                trainingData.add( datum );
            }
        }
        logger.info( "Training with " + trainingData.size() + " variants after standard deviation thresholding." );
        if( trainingData.size() < VRAC.MIN_NUM_BAD_VARIANTS ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
        }
        return trainingData;
    }

    public ExpandingArrayList<VariantDatum> selectWorstVariants( double bottomPercentage, final int minimumNumber ) {
        // The return value is the list of training variants
        final ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();

        // First add to the training list all sites overlapping any bad sites training tracks
        for( final VariantDatum datum : data ) {
            if( datum.atAntiTrainingSite && !datum.failingSTDThreshold && !Double.isInfinite(datum.lod) ) {
                trainingData.add( datum );
            }
        }
        final int numBadSitesAdded = trainingData.size();
        logger.info( "Found " + numBadSitesAdded + " variants overlapping bad sites training tracks." );

        // Next sort the variants by the LOD coming from the positive model and add to the list the bottom X percent of variants
        Collections.sort( data );
        final int numToAdd = Math.max( minimumNumber - trainingData.size(), Math.round((float)bottomPercentage * data.size()) );
        if( numToAdd > data.size() ) {
            throw new UserException.BadInput( "Error during negative model training. Minimum number of variants to use in training is larger than the whole call set. One can attempt to lower the --minNumBadVariants arugment but this is unsafe." );
        } else if( numToAdd == minimumNumber - trainingData.size() ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
            bottomPercentage = ((float) numToAdd) / ((float) data.size());
        }
        int index = 0, numAdded = 0;
        while( numAdded < numToAdd ) {
            final VariantDatum datum = data.get(index++);
            if( !datum.atAntiTrainingSite && !datum.failingSTDThreshold && !Double.isInfinite(datum.lod) ) {
                datum.atAntiTrainingSite = true;
                trainingData.add( datum );
                numAdded++;
            }
        }
        logger.info( "Additionally training with worst " + String.format("%.3f", (float) bottomPercentage * 100.0f) + "% of passing data --> " + (trainingData.size() - numBadSitesAdded) + " variants with LOD <= " + String.format("%.4f", data.get(index).lod) + "." );
        return trainingData;
    }

    public ExpandingArrayList<VariantDatum> getRandomDataForPlotting( int numToAdd ) {
        numToAdd = Math.min(numToAdd, data.size());
        final ExpandingArrayList<VariantDatum> returnData = new ExpandingArrayList<VariantDatum>();
        for( int iii = 0; iii < numToAdd; iii++) {
            final VariantDatum datum = data.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(data.size()));
            if( !datum.failingSTDThreshold ) {
                returnData.add(datum);
            }
        }

        // Add an extra 5% of points from bad training set, since that set is small but interesting
        for( int iii = 0; iii < Math.floor(0.05*numToAdd); iii++) {
            final VariantDatum datum = data.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(data.size()));
            if( datum.atAntiTrainingSite && !datum.failingSTDThreshold ) { returnData.add(datum); }
            else { iii--; }
        }

        return returnData;
    }

    private double mean( final int index ) {
        double sum = 0.0;
        int numNonNull = 0;
        for( final VariantDatum datum : data ) {
            if( datum.atTrainingSite && !datum.isNull[index] ) { sum += datum.annotations[index]; numNonNull++; }
        }
        return sum / ((double) numNonNull);
    }

    private double standardDeviation( final double mean, final int index ) {
        double sum = 0.0;
        int numNonNull = 0;
        for( final VariantDatum datum : data ) {
            if( datum.atTrainingSite && !datum.isNull[index] ) { sum += ((datum.annotations[index] - mean)*(datum.annotations[index] - mean)); numNonNull++; }
        }
        return Math.sqrt( sum / ((double) numNonNull) );
    }

    public void decodeAnnotations( final VariantDatum datum, final VariantContext vc, final boolean jitter ) {
        final double[] annotations = new double[annotationKeys.size()];
        final boolean[] isNull = new boolean[annotationKeys.size()];
        int iii = 0;
        for( final String key : annotationKeys ) {
            isNull[iii] = false;
            annotations[iii] = decodeAnnotation( key, vc, jitter );
            if( Double.isNaN(annotations[iii]) ) { isNull[iii] = true; }
            iii++;
        }
        datum.annotations = annotations;
        datum.isNull = isNull;
    }

    private static double decodeAnnotation( final String annotationKey, final VariantContext vc, final boolean jitter ) {
        double value;

        try {
            value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            if( Double.isInfinite(value) ) { value = Double.NaN; }
            if( jitter && annotationKey.equalsIgnoreCase("HRUN") ) { // Integer valued annotations must be jittered a bit to work in this GMM
                  value += -0.25 + 0.5 * GenomeAnalysisEngine.getRandomGenerator().nextDouble();
            }

            if (vc.isIndel() && annotationKey.equalsIgnoreCase("QD")) {
            // normalize QD by event length for indel case
                int eventLength = Math.abs(vc.getAlternateAllele(0).getBaseString().length() - vc.getReference().getBaseString().length()); // ignore multi-allelic complication here for now
                if (eventLength > 0) { // sanity check
                    value /= (double)eventLength;
                }
            }

            if( jitter && annotationKey.equalsIgnoreCase("HaplotypeScore") && MathUtils.compareDoubles(value, 0.0, 0.0001) == 0 ) { value = -0.2 + 0.4*GenomeAnalysisEngine.getRandomGenerator().nextDouble(); }
            if( jitter && annotationKey.equalsIgnoreCase("FS") && MathUtils.compareDoubles(value, 0.0, 0.001) == 0 ) { value = -0.2 + 0.4*GenomeAnalysisEngine.getRandomGenerator().nextDouble(); }
        } catch( Exception e ) {
            value = Double.NaN; // The VQSR works with missing data by marginalizing over the missing dimension when evaluating the Gaussian mixture model
        }

        return value;
    }

    public void parseTrainingSets( final RefMetaDataTracker tracker, final GenomeLoc genomeLoc, final VariantContext evalVC, final VariantDatum datum, final boolean TRUST_ALL_POLYMORPHIC ) {
        datum.isKnown = false;
        datum.atTruthSite = false;
        datum.atTrainingSite = false;
        datum.atAntiTrainingSite = false;
        datum.prior = 2.0;

        for( final TrainingSet trainingSet : trainingSets ) {
            for( final VariantContext trainVC : tracker.getValues(trainingSet.rodBinding, genomeLoc) ) {
                if( isValidVariant( evalVC, trainVC, TRUST_ALL_POLYMORPHIC ) ) {
                    datum.isKnown = datum.isKnown || trainingSet.isKnown;
                    datum.atTruthSite = datum.atTruthSite || trainingSet.isTruth;
                    datum.atTrainingSite = datum.atTrainingSite || trainingSet.isTraining;
                    datum.prior = Math.max( datum.prior, trainingSet.prior );
                    datum.consensusCount += ( trainingSet.isConsensus ? 1 : 0 );
                }
                if( trainVC != null ) {
                    datum.atAntiTrainingSite = datum.atAntiTrainingSite || trainingSet.isAntiTraining;
                }
            }
        }
    }

    private boolean isValidVariant( final VariantContext evalVC, final VariantContext trainVC, final boolean TRUST_ALL_POLYMORPHIC) {
        return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() &&
                        ((evalVC.isSNP() && trainVC.isSNP()) || ((evalVC.isIndel()||evalVC.isMixed()) && (trainVC.isIndel()||trainVC.isMixed()))) &&
                        (TRUST_ALL_POLYMORPHIC || !trainVC.hasGenotypes() || trainVC.isPolymorphic());
    }

    public void writeOutRecalibrationTable( final PrintStream RECAL_FILE ) {
        for( final VariantDatum datum : data ) {
            RECAL_FILE.println(String.format("%s,%d,%d,%.4f,%s",
                    datum.contig, datum.start, datum.stop, datum.lod,
                    (datum.worstAnnotation != -1 ? annotationKeys.get(datum.worstAnnotation) : "NULL")));
        }
    }
}

package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * Date: Nov 6, 2009
 *
 * This helper class holds the data HashMap as well as submaps that represent the marginal distributions collapsed over all needed dimensions.
 */

public class RecalDataManager {
    public NHashMap<RecalDatum> data; // The full dataset
    private NHashMap<RecalDatum> dataCollapsedReadGroup; // Table where everything except read group has been collapsed
    private NHashMap<RecalDatum> dataCollapsedQualityScore; // Table where everything except read group and quality score has been collapsed
    private ArrayList<NHashMap<RecalDatum>> dataCollapsedByCovariate; // Tables where everything except read group, quality score, and given covariate has been collapsed
    private NHashMap<Double> dataSumExpectedErrors; // Table used to calculate the overall aggregate reported quality score in which everything except read group is collapsed

    private NHashMap<Double> dataCollapsedReadGroupDouble; // Table of empirical qualities where everything except read group has been collapsed
    private NHashMap<Double> dataCollapsedQualityScoreDouble; // Table of empirical qualities where everything except read group and quality score has been collapsed
    private ArrayList<NHashMap<Double>> dataCollapsedByCovariateDouble; // Table of empirical qualities where everything except read group, quality score, and given covariate has been collapsed
    public NHashMap<Double> aggregateReportedQuality; // Table of the overall aggregate reported quality score in which everything except read group is collapsed

    public final static String ORIGINAL_QUAL_ATTRIBUTE_TAG = "OQ"; // The tag that holds the original quality scores
    public final static String COLOR_SPACE_QUAL_ATTRIBUTE_TAG = "CQ"; // The tag that holds the color space quality scores for SOLID bams
    private static boolean warnUserNullReadGroup = false;

    RecalDataManager() {
    	data = new NHashMap<RecalDatum>();
    }

    RecalDataManager( final int estimatedCapacity, final boolean createCollapsedTables, final int numCovariates ) {
    	if( createCollapsedTables ) { // Initialize all the collapsed tables, only used by TableRecalibrationWalker
            dataCollapsedReadGroup = new NHashMap<RecalDatum>();
            dataCollapsedQualityScore = new NHashMap<RecalDatum>();
            dataCollapsedByCovariate = new ArrayList<NHashMap<RecalDatum>>();
            for( int iii = 0; iii < numCovariates - 2; iii++ ) { // readGroup and QualityScore aren't counted here, their tables are separate
                dataCollapsedByCovariate.add( new NHashMap<RecalDatum>() );
            }
            dataSumExpectedErrors = new NHashMap<Double>();
            aggregateReportedQuality = new NHashMap<Double>();
        } else {
            data = new NHashMap<RecalDatum>( estimatedCapacity, 0.8f);
        }            
    }
    
    RecalDataManager( final int estimatedCapacity ) {
        data = new NHashMap<RecalDatum>( estimatedCapacity, 0.8f ); // Second arg is the 'loading factor',
                                                                    //   a number to monkey around with when optimizing performace of the HashMap
    }

    /**
     * Add the given mapping to all of the collapsed hash tables
     * @param key The list of comparables that is the key for this mapping
     * @param fullDatum The RecalDatum which is the data for this mapping
     */
    public final void addToAllTables( final List<? extends Comparable> key, final RecalDatum fullDatum ) {

        // The full dataset isn't actually ever used for anything because of the sequential calculation so no need to keep the full data HashMap around
        //data.put(key, thisDatum); // add the mapping to the main table

        // Create dataCollapsedReadGroup, the table where everything except read group has been collapsed
        //BUGBUG: Bases with Reported Quality less than -pQ argument shouldn't be included in this collapsed table
        ArrayList<Comparable> newKey = new ArrayList<Comparable>();
        newKey.add( key.get(0) ); // Make a new key with just the read group
        RecalDatum collapsedDatum = dataCollapsedReadGroup.get( newKey );
        if( collapsedDatum == null ) {
            dataCollapsedReadGroup.put( newKey, new RecalDatum(fullDatum) );
        } else {
            collapsedDatum.increment(fullDatum);
        }

        // Create dataSumExpectedErrors, the table used to calculate the overall aggregate quality score in which everything except read group is collapsed
        newKey = new ArrayList<Comparable>();
        newKey.add( key.get(0) ); // Make a new key with just the read group
        Double sumExpectedErrors = dataSumExpectedErrors.get( newKey );
        if( sumExpectedErrors == null ) {
            dataSumExpectedErrors.put( newKey, 0.0 );
        } else {
            dataSumExpectedErrors.remove( newKey );
            sumExpectedErrors += QualityUtils.qualToErrorProb(Byte.parseByte(key.get(1).toString())) * fullDatum.getNumObservations();
            dataSumExpectedErrors.put( newKey, sumExpectedErrors );
        }

        newKey = new ArrayList<Comparable>();
        // Create dataCollapsedQuality, the table where everything except read group and quality score has been collapsed
        newKey.add( key.get(0) ); // Make a new key with the read group ...
        newKey.add( key.get(1) ); //                                    and quality score
        collapsedDatum = dataCollapsedQualityScore.get( newKey );
        if( collapsedDatum == null ) {
            dataCollapsedQualityScore.put( newKey, new RecalDatum(fullDatum) );
        } else {
            collapsedDatum.increment(fullDatum);
        }

        // Create dataCollapsedByCovariate's, the tables where everything except read group, quality score, and given covariate has been collapsed
        for( int iii = 0; iii < dataCollapsedByCovariate.size(); iii++ ) {
            newKey = new ArrayList<Comparable>();
            newKey.add( key.get(0) ); // Make a new key with the read group ...
            newKey.add( key.get(1) ); //                                    and quality score ...
            newKey.add( key.get(iii + 2) ); //                                                    and the given covariate
            collapsedDatum = dataCollapsedByCovariate.get(iii).get( newKey );
            if( collapsedDatum == null ) {
                dataCollapsedByCovariate.get(iii).put( newKey, new RecalDatum(fullDatum) );
            } else {
                collapsedDatum.increment(fullDatum);
            }
        }
    }

    /**
     * Loop over all the collapsed tables and turn the recalDatums found there into an empricial quality score
     *   that will be used in the sequential calculation in TableRecalibrationWalker
     * @param numCovariates The number of covariates you have determines the number of tables to create
     * @param smoothing The smoothing paramter that goes into empirical quality score calculation
     */
    public final void generateEmpiricalQualities( final int numCovariates, final int smoothing ) {

        dataCollapsedReadGroupDouble = new NHashMap<Double>();
        dataCollapsedQualityScoreDouble = new NHashMap<Double>();
        dataCollapsedByCovariateDouble = new ArrayList<NHashMap<Double>>();
        for( int iii = 0; iii < numCovariates - 2; iii++ ) { // readGroup and QualityScore aren't counted here, their tables are separate
            dataCollapsedByCovariateDouble.add( new NHashMap<Double>() );
        }

        // Hash the empirical quality scores so we don't have to call Math.log at every base for every read
        // Looping over the entrySet is really expensive but worth it
        for( Map.Entry<List<? extends Comparable>,RecalDatum> entry : dataCollapsedReadGroup.entrySet() ) {
            dataCollapsedReadGroupDouble.put( entry.getKey(), entry.getValue().empiricalQualDouble( smoothing ) );
            aggregateReportedQuality.put( entry.getKey(), // sum of the expected errors divided by number of observations turned into a Q score
                    QualityUtils.phredScaleErrorRate( dataSumExpectedErrors.get(entry.getKey()) / ((double)entry.getValue().getNumObservations()) ) );
        }
        for( Map.Entry<List<? extends Comparable>,RecalDatum> entry : dataCollapsedQualityScore.entrySet() ) {
            dataCollapsedQualityScoreDouble.put( entry.getKey(), entry.getValue().empiricalQualDouble( smoothing ) );
        }
        for( int iii = 0; iii < numCovariates - 2; iii++ ) {
            for( Map.Entry<List<? extends Comparable>,RecalDatum> entry : dataCollapsedByCovariate.get(iii).entrySet() ) {
                dataCollapsedByCovariateDouble.get(iii).put( entry.getKey(), entry.getValue().empiricalQualDouble( smoothing ) );
            }
        }

        dataSumExpectedErrors.clear();
        dataCollapsedReadGroup.clear();
        dataCollapsedQualityScore.clear();
        for( int iii = 0; iii < numCovariates - 2; iii++ ) {
            dataCollapsedByCovariate.get(iii).clear();
        }
        dataCollapsedByCovariate.clear();
        dataSumExpectedErrors = null; // Will never need this table again
        dataCollapsedReadGroup = null; // Will never need this table again
        dataCollapsedQualityScore = null; // Will never need this table again
        dataCollapsedByCovariate = null; // Will never need this table again
        if( data != null ) {
            data.clear();
            data = null; // Will never need this table again
        }
    }

    /**
     * Get the appropriate collapsed table out of the set of all the tables held by this Object
     * @param covariate Which covariate indexes the desired collapsed HashMap
     * @return The desired collapsed HashMap
     */
    public final NHashMap<RecalDatum> getCollapsedTable( final int covariate ) {
        if( covariate == 0) {
            return dataCollapsedReadGroup; // Table where everything except read group has been collapsed
        } else if( covariate == 1 ) {
            return dataCollapsedQualityScore; // Table where everything except read group and quality score has been collapsed
        } else {
            return dataCollapsedByCovariate.get( covariate - 2 ); // Table where everything except read group, quality score, and given covariate has been collapsed
        }
    }

    /**
     * Get the appropriate collapsed table of emprical quality out of the set of all the tables held by this Object
     * @param covariate Which covariate indexes the desired collapsed NHashMap<Double>
     * @return The desired collapsed NHashMap<Double>
     */
    public final NHashMap<Double> getCollapsedDoubleTable( final int covariate ) {
        if( covariate == 0) {
            return dataCollapsedReadGroupDouble; // Table of empirical qualities where everything except read group has been collapsed
        } else if( covariate == 1 ) {
            return dataCollapsedQualityScoreDouble; // Table of empirical qualities where everything except read group and quality score has been collapsed
        } else {
            return dataCollapsedByCovariateDouble.get( covariate - 2 ); // Table of empirical qualities where everything except read group, quality score, and given covariate has been collapsed
        }
    }

    /**
     * Section of code shared between the two recalibration walkers which uses the command line arguments to adjust attributes of the read such as quals or platform string
     * @param read The read to adjust
     * @param RAC The list of shared command line arguments
     */
    public static void parseSAMRecord( final SAMRecord read, final RecalibrationArgumentCollection RAC) {

        // Check if we need to use the original quality scores instead
        if ( RAC.USE_ORIGINAL_QUALS && read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG) != null ) {
            Object obj = read.getAttribute(RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG);
            if ( obj instanceof String )
                read.setBaseQualities( QualityUtils.fastqToPhred((String)obj) );
            else {
                throw new RuntimeException(String.format("Value encoded by %s in %s isn't a string!", RecalDataManager.ORIGINAL_QUAL_ATTRIBUTE_TAG, read.getReadName()));
            }
        }

        SAMReadGroupRecord readGroup = read.getReadGroup();

        // If there are no read groups we have to default to something, and that something could be specified by the user using command line arguments
        if( readGroup == null ) {
            if( !warnUserNullReadGroup && RAC.FORCE_READ_GROUP == null ) {
                Utils.warnUser("The input .bam file contains reads with no read group. " +
                                "Defaulting to read group ID = " + RAC.DEFAULT_READ_GROUP + " and platform = " + RAC.DEFAULT_PLATFORM + ". " +
                                "First observed at read with name = " + read.getReadName() );
                warnUserNullReadGroup = true;
            }
            // There is no readGroup so defaulting to these values
            readGroup = new SAMReadGroupRecord( RAC.DEFAULT_READ_GROUP );
            readGroup.setPlatform( RAC.DEFAULT_PLATFORM );
            ((GATKSAMRecord)read).setReadGroup( readGroup );
        }

        if( RAC.FORCE_READ_GROUP != null && !read.getReadGroup().getReadGroupId().equals(RAC.FORCE_READ_GROUP) ) { // Collapse all the read groups into a single common String provided by the user
            String oldPlatform = readGroup.getPlatform();
            readGroup = new SAMReadGroupRecord( RAC.FORCE_READ_GROUP );
            readGroup.setPlatform( oldPlatform );
            ((GATKSAMRecord)read).setReadGroup( readGroup );
        }

        if( RAC.FORCE_PLATFORM != null && !read.getReadGroup().getPlatform().equals(RAC.FORCE_PLATFORM)) {
            read.getReadGroup().setPlatform( RAC.FORCE_PLATFORM );
        }
    }
}

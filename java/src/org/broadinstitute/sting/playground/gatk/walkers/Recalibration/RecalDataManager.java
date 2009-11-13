package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.*;

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
 */

public class RecalDataManager {
    public NHashMap<RecalDatum> data; // the full dataset
    private NHashMap<RecalDatum> dataCollapsedReadGroup; // table where everything except read group has been collapsed
    private NHashMap<RecalDatum> dataCollapsedQualityScore; // table where everything except read group and quality score has been collapsed
    private ArrayList<NHashMap<RecalDatum>> dataCollapsedByCovariate; // tables where everything except read group, quality score, and given covariate has been collapsed
    private boolean collapsedTablesCreated;
    public NHashMap<Double> dataSumExpectedErrors;

    public final static String ORIGINAL_QUAL_ATTRIBUTE_TAG = "OQ";

    RecalDataManager() {
    	data = new NHashMap<RecalDatum>();
    	collapsedTablesCreated = false;
    }
    
    RecalDataManager( int estimatedCapacity ) {
        data = new NHashMap<RecalDatum>( estimatedCapacity, 0.95f ); // second arg is the 'loading factor',
                                                                     //   a number to monkey around with when optimizing performace of the HashMap
        collapsedTablesCreated = false;
    }

    // BUGBUG: A lot going on in this method, doing a lot of pre-calculations for use in the sequential mode calculation later in TableRecalibrationWalker
    public final void createCollapsedTables( final int numCovariates ) {
        dataCollapsedReadGroup = new NHashMap<RecalDatum>();
        dataCollapsedQualityScore = new NHashMap<RecalDatum>();
        dataCollapsedByCovariate = new ArrayList<NHashMap<RecalDatum>>();
        for( int iii = 0; iii < numCovariates - 2; iii++ ) { // readGroup and QualityScore aren't counted here, their tables are separate
            dataCollapsedByCovariate.add( new NHashMap<RecalDatum>() );
        }
        dataSumExpectedErrors = new NHashMap<Double>();

        // preallocate for use in for loops below
        RecalDatum thisDatum;
        RecalDatum collapsedDatum;
        List<? extends Comparable> key;
        ArrayList<Comparable> newKey;
        Double sumExpectedErrors;

        // for every data point in the map
        for( Map.Entry<List<? extends Comparable>,RecalDatum> entry : data.entrySet() ) {
            thisDatum = entry.getValue();
            key = entry.getKey();
            
            // create dataCollapsedReadGroup, the table where everything except read group has been collapsed
            newKey = new ArrayList<Comparable>();
            newKey.add( key.get(0) ); // make a new key with just the read group
            collapsedDatum = dataCollapsedReadGroup.get( newKey );
            if( collapsedDatum == null ) {
                dataCollapsedReadGroup.put( newKey, new RecalDatum( thisDatum ) );
            } else {
                collapsedDatum.increment( thisDatum );
            }
            
            newKey = new ArrayList<Comparable>();
            newKey.add( key.get(0) ); // make a new key with just the read group
            sumExpectedErrors = dataSumExpectedErrors.get( newKey );
            if( sumExpectedErrors == null ) {
                dataSumExpectedErrors.put( newKey, 0.0 );
            } else {
                dataSumExpectedErrors.remove( newKey );
                sumExpectedErrors += QualityUtils.qualToErrorProb(Byte.parseByte(key.get(1).toString())) * thisDatum.getNumObservations();
                dataSumExpectedErrors.put( newKey, sumExpectedErrors );
            }

            newKey = new ArrayList<Comparable>();
            // create dataCollapsedQuality, the table where everything except read group and quality score has been collapsed
            newKey.add( key.get(0) ); // make a new key with the read group ...
            newKey.add( key.get(1) ); //                                    and quality score
            collapsedDatum = dataCollapsedQualityScore.get( newKey );
            if( collapsedDatum == null ) {
            	dataCollapsedQualityScore.put( newKey, new RecalDatum( thisDatum ) );
            } else {
                collapsedDatum.increment( thisDatum );
            }

            // create dataCollapsedByCovariate's, the tables where everything except read group, quality score, and given covariate has been collapsed
            for( int iii = 0; iii < numCovariates - 2; iii++ ) { // readGroup and QualityScore aren't counted
                newKey = new ArrayList<Comparable>();
                newKey.add( key.get(0) ); // make a new key with the read group ...
                newKey.add( key.get(1) ); //                                    and quality score ...
                newKey.add( key.get(iii + 2) ); //                                                    and the given covariate
                collapsedDatum = dataCollapsedByCovariate.get(iii).get( newKey );
                if( collapsedDatum == null ) {
                    dataCollapsedByCovariate.get(iii).put( newKey, new RecalDatum( thisDatum ) );
                } else {
                    collapsedDatum.increment( thisDatum );
                }
            }
        }

        collapsedTablesCreated = true;
    }

    public final NHashMap<RecalDatum> getCollapsedTable( final int covariate ) {
        if( !collapsedTablesCreated ) {
            throw new StingException("Trying to get collapsed tables before they have been populated. Null pointers abound.");
        }

        if( covariate == 0) {
            return dataCollapsedReadGroup; // table where everything except read group has been collapsed
        } else if( covariate == 1 ) {
            return dataCollapsedQualityScore; // table where everything except read group and quality score has been collapsed
        } else {
            return dataCollapsedByCovariate.get( covariate - 2 ); // table where everything except read group, quality score, and given covariate has been collapsed
        }
    }
}

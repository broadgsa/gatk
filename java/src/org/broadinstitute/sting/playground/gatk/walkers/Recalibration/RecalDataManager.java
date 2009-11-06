package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 6, 2009
 */

public class RecalDataManager {
    public NHashMap<RecalDatum> data; // the full dataset
    public NHashMap<RecalDatum> dataCollapsedReadGroup; // table where everything except read group has been collapsed
    public NHashMap<RecalDatum> dataCollapsedQualityScore; // table where everything except read group and quality score has been collapsed
    public ArrayList<NHashMap<RecalDatum>> dataCollapsedByCovariate; // tables where everything except read group, quality score, and given covariate has been collapsed
    public boolean collapsedTablesCreated;
    public NHashMap<Double> dataSumExpectedErrors;

    RecalDataManager() {
        data = new NHashMap<RecalDatum>();
        collapsedTablesCreated = false;
    }

    // BUGBUG: A lot going on in this method, doing a lot of pre-calculations for use in the sequential mode calculation later
    public void createCollapsedTables( int numCovariates ) {
        dataCollapsedReadGroup = new NHashMap<RecalDatum>();
        dataCollapsedQualityScore = new NHashMap<RecalDatum>();
        dataCollapsedByCovariate = new ArrayList<NHashMap<RecalDatum>>();
        for( int iii = 0; iii < numCovariates - 2; iii++ ) { // readGroup and QualityScore aren't counted
            dataCollapsedByCovariate.add( new NHashMap<RecalDatum>() );
        }
        dataSumExpectedErrors = new NHashMap<Double>();

        // preallocate for use in for loops below
        RecalDatum thisDatum;
        RecalDatum collapsedDatum;
        List<? extends Comparable<?>> key;
        ArrayList<Comparable<?>> newKey;
        Double sumExpectedErrors;

        // for every data point in the map
        for( Map.Entry<List<? extends Comparable<?>>,RecalDatum> entry : data.entrySet() ) {
            thisDatum = entry.getValue();
            key = entry.getKey();
            
            // create dataCollapsedReadGroup, the table where everything except read group has been collapsed
            newKey = new ArrayList<Comparable<?>>();
            newKey.add( key.get(0) ); // make a new key with just the read group
            collapsedDatum = dataCollapsedReadGroup.get( newKey );
            if( collapsedDatum == null ) {
                dataCollapsedReadGroup.put( newKey, new RecalDatum( thisDatum ) );
                //System.out.println("Added key: " + newKey + " to the dataCollapsedReadGroup");
            } else {
                collapsedDatum.increment( thisDatum );
            }

            newKey = new ArrayList<Comparable<?>>();
            newKey.add( key.get(0) ); // make a new key with just the read group
            sumExpectedErrors = dataSumExpectedErrors.get( newKey );
            if( sumExpectedErrors == null ) {
                dataSumExpectedErrors.put( newKey, 0.0 );
            } else {
                //System.out.println("updated += " + QualityUtils.qualToErrorProb(Byte.parseByte(key.get(1).toString())) * thisDatum.getNumObservations());
                dataSumExpectedErrors.remove( newKey );
                sumExpectedErrors += QualityUtils.qualToErrorProb(Byte.parseByte(key.get(1).toString())) * thisDatum.getNumObservations();
                dataSumExpectedErrors.put( newKey, sumExpectedErrors );
            }

            newKey = new ArrayList<Comparable<?>>();
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
                newKey = new ArrayList<Comparable<?>>();
                newKey.add( key.get(0) ); // make a new key with the read group ...
                newKey.add( key.get(1) ); //                                    and quality score ...
                newKey.add( key.get(iii) ); //                                                    and the given covariate
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

    public NHashMap<RecalDatum> getCollapsedTable( int covariate ) {
        if( !collapsedTablesCreated ) {
            throw new StingException("Trying to get collapsed tables before they have been populated.");
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

package org.broadinstitute.sting.analyzecovariates;

import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum;
import org.broadinstitute.sting.gatk.walkers.recalibration.NHashMap;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Dec 1, 2009
 *
 * The difference between this AnalysisDataManager and the RecalDataManager used by the Recalibration walkers is that here the collapsed data tables are indexed
 *  by only read group and the given covariate, while in the recalibrator the collapsed tables are indexed by read group, reported quality, and the given covariate.
 */

public class AnalysisDataManager {
    
    private NHashMap<RecalDatum> dataCollapsedReadGroup; // Table where everything except read group has been collapsed
    private ArrayList<NHashMap<RecalDatum>> dataCollapsedByCovariate; // Tables where everything except read group and given covariate has been collapsed

    AnalysisDataManager() {
    }

    AnalysisDataManager( final int numCovariates ) {
        dataCollapsedReadGroup = new NHashMap<RecalDatum>();
        dataCollapsedByCovariate = new ArrayList<NHashMap<RecalDatum>>();
        for( int iii = 0; iii < numCovariates - 1; iii++ ) { // readGroup isn't counted here, its table is separate
            dataCollapsedByCovariate.add( new NHashMap<RecalDatum>() );
        }
    }

    /**
     * Add the given mapping to all of the collapsed hash tables
     * @param key The list of comparables that is the key for this mapping
     * @param fullDatum The RecalDatum which is the data for this mapping
     */
    public final void addToAllTables( final List<? extends Comparable> key, final RecalDatum fullDatum, final int IGNORE_QSCORES_LESS_THAN ) {

        int qscore = Integer.parseInt( key.get(1).toString() );
        ArrayList<Comparable> newKey;
        RecalDatum collapsedDatum;

        if( !(qscore < IGNORE_QSCORES_LESS_THAN) ) {
            // Create dataCollapsedReadGroup, the table where everything except read group has been collapsed
            newKey = new ArrayList<Comparable>();
            newKey.add( key.get(0) ); // Make a new key with just the read group
            collapsedDatum = dataCollapsedReadGroup.get( newKey );
            if( collapsedDatum == null ) {
                dataCollapsedReadGroup.put( newKey, new RecalDatum(fullDatum) );
            } else {
                collapsedDatum.combine( fullDatum ); // using combine instead of increment in order to calculate overall aggregateQReported
            }
        }

        // Create dataCollapsedByCovariate's, the tables where everything except read group, quality score, and given covariate has been collapsed
        for( int iii = 0; iii < dataCollapsedByCovariate.size(); iii++ ) {
            if( iii == 0 || !(qscore < IGNORE_QSCORES_LESS_THAN) ) { // use all data for the plot versus reported quality, but not for the other plots versus cycle and etc.
                newKey = new ArrayList<Comparable>();
                newKey.add( key.get(0) ); // Make a new key with the read group ...
                newKey.add( key.get(iii + 1) ); //                                                    and the given covariate
                collapsedDatum = dataCollapsedByCovariate.get(iii).get( newKey );
                if( collapsedDatum == null ) {
                    dataCollapsedByCovariate.get(iii).put( newKey, new RecalDatum(fullDatum) );
                } else {
                    collapsedDatum.combine( fullDatum );
                }
            }
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
        } else {
            return dataCollapsedByCovariate.get( covariate - 1 ); // Table where everything except read group, quality score, and given covariate has been collapsed
        }
    }

}
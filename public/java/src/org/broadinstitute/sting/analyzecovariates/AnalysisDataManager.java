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

package org.broadinstitute.sting.analyzecovariates;

import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.collections.NestedHashMap;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Dec 1, 2009
 *
 * The difference between this AnalysisDataManager and the RecalDataManager used by the Recalibration walkers is that here the collapsed data tables are indexed
 *  by only read group and the given covariate, while in the recalibrator the collapsed tables are indexed by read group, reported quality, and the given covariate.
 */

public class AnalysisDataManager {
    
    private NestedHashMap dataCollapsedReadGroup; // Table where everything except read group has been collapsed
    private ArrayList<NestedHashMap> dataCollapsedByCovariate; // Tables where everything except read group and given covariate has been collapsed

    AnalysisDataManager() {
    }

    AnalysisDataManager( final int numCovariates ) {
        dataCollapsedReadGroup = new NestedHashMap();
        dataCollapsedByCovariate = new ArrayList<NestedHashMap>();
        for( int iii = 0; iii < numCovariates - 1; iii++ ) { // readGroup isn't counted here, its table is separate
            dataCollapsedByCovariate.add( new NestedHashMap() );
        }
    }

    /**
     * Add the given mapping to all of the collapsed hash tables
     * @param key The list of comparables that is the key for this mapping
     * @param fullDatum The RecalDatum which is the data for this mapping
     * @param IGNORE_QSCORES_LESS_THAN The threshold in report quality for adding to the aggregate collapsed table
     */
    public final void addToAllTables( final Object[] key, final RecalDatum fullDatum, final int IGNORE_QSCORES_LESS_THAN ) {

        int qscore = Integer.parseInt( key[1].toString() );
        RecalDatum collapsedDatum;
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] covariateCollapsedKey = new Object[2];

        if( !(qscore < IGNORE_QSCORES_LESS_THAN) ) {
            // Create dataCollapsedReadGroup, the table where everything except read group has been collapsed
            readGroupCollapsedKey[0] = key[0]; // Make a new key with just the read group
            collapsedDatum = (RecalDatum)dataCollapsedReadGroup.get( readGroupCollapsedKey );
            if( collapsedDatum == null ) {
                dataCollapsedReadGroup.put( new RecalDatum(fullDatum), readGroupCollapsedKey );
            } else {
                collapsedDatum.combine( fullDatum ); // using combine instead of increment in order to calculate overall aggregateQReported
            }
        }

        // Create dataCollapsedByCovariate's, the tables where everything except read group and given covariate has been collapsed
        for( int iii = 0; iii < dataCollapsedByCovariate.size(); iii++ ) {
            if( iii == 0 || !(qscore < IGNORE_QSCORES_LESS_THAN) ) { // use all data for the plot versus reported quality, but not for the other plots versus cycle and etc.
                covariateCollapsedKey[0] = key[0]; // Make a new key with the read group ...
                Object theCovariateElement = key[iii + 1]; //           and the given covariate
                if( theCovariateElement != null ) {
                    covariateCollapsedKey[1] = theCovariateElement;
                    collapsedDatum = (RecalDatum)dataCollapsedByCovariate.get(iii).get( covariateCollapsedKey );
                    if( collapsedDatum == null ) {
                        dataCollapsedByCovariate.get(iii).put( new RecalDatum(fullDatum), covariateCollapsedKey );
                    } else {
                        collapsedDatum.combine( fullDatum );
                    }
                }
            }
        }
    }

    /**
     * Get the appropriate collapsed table out of the set of all the tables held by this Object
     * @param covariate Which covariate indexes the desired collapsed HashMap
     * @return The desired collapsed HashMap
     */
    public final NestedHashMap getCollapsedTable( final int covariate ) {
        if( covariate == 0) {
            return dataCollapsedReadGroup; // Table where everything except read group has been collapsed
        } else {
            return dataCollapsedByCovariate.get( covariate - 1 ); // Table where everything except read group, quality score, and given covariate has been collapsed
        }
    }

}
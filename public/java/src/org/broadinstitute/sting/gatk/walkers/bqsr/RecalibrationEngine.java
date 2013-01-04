/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;
import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;

public class RecalibrationEngine {
    final protected Covariate[] covariates;
    final private int numReadGroups;
    final private PrintStream maybeLogStream;
    final private boolean lowMemoryMode;

    /**
     * Has finalizeData() been called?
     */
    private boolean finalized = false;

    /**
     * The final (merged, etc) recalibration tables, suitable for downstream analysis.
     */
    private RecalibrationTables finalRecalibrationTables = null;

    private final List<RecalibrationTables> recalibrationTablesList = new LinkedList<RecalibrationTables>();

    private final ThreadLocal<RecalibrationTables> threadLocalTables = new ThreadLocal<RecalibrationTables>() {
        private synchronized RecalibrationTables makeAndCaptureTable() {
            final RecalibrationTables newTable = new RecalibrationTables(covariates, numReadGroups, maybeLogStream);
            recalibrationTablesList.add(newTable);
            return newTable;
        }

        @Override
        protected synchronized RecalibrationTables initialValue() {
            if ( lowMemoryMode ) {
                return recalibrationTablesList.isEmpty() ? makeAndCaptureTable() : recalibrationTablesList.get(0);
            } else {
                return makeAndCaptureTable();
            }
        }
    };

    /**
     * Get a recalibration table suitable for updating the underlying RecalDatums
     *
     * May return a thread-local version, or a single version, depending on the initialization
     * arguments of this instance.
     *
     * @return updated tables
     */
    protected RecalibrationTables getUpdatableRecalibrationTables() {
        return threadLocalTables.get();
    }

    /**
     * Initialize the recalibration engine
     *
     * Called once before any calls to updateDataForRead are made.  The engine should prepare itself
     * to handle any number of updateDataForRead calls containing ReadRecalibrationInfo containing
     * keys for each of the covariates provided.
     *
     * The engine should collect match and mismatch data into the recalibrationTables data.
     *
     * @param covariates an array of the covariates we'll be using in this engine, order matters
     * @param numReadGroups the number of read groups we should use for the recalibration tables
     * @param maybeLogStream an optional print stream for logging calls to the nestedhashmap in the recalibration tables
     */
    public RecalibrationEngine(final Covariate[] covariates, final int numReadGroups, final PrintStream maybeLogStream, final boolean enableLowMemoryMode) {
        if ( covariates == null ) throw new IllegalArgumentException("Covariates cannot be null");
        if ( numReadGroups < 1 ) throw new IllegalArgumentException("numReadGroups must be >= 1 but got " + numReadGroups);

        this.covariates = covariates.clone();
        this.numReadGroups = numReadGroups;
        this.maybeLogStream = maybeLogStream;
        this.lowMemoryMode = enableLowMemoryMode;
    }

    /**
     * Update the recalibration statistics using the information in recalInfo
     * @param recalInfo data structure holding information about the recalibration values for a single read
     */
    @Requires("recalInfo != null")
    public void updateDataForRead( final ReadRecalibrationInfo recalInfo ) {
        final GATKSAMRecord read = recalInfo.getRead();
        final ReadCovariates readCovariates = recalInfo.getCovariatesValues();
        final RecalibrationTables tables = getUpdatableRecalibrationTables();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = tables.getQualityScoreTable();

        for( int offset = 0; offset < read.getReadBases().length; offset++ ) {
            if( ! recalInfo.skip(offset) ) {

                for (final EventType eventType : EventType.values()) {
                    final int[] keys = readCovariates.getKeySet(offset, eventType);
                    final int eventIndex = eventType.ordinal();
                    final byte qual = recalInfo.getQual(eventType, offset);
                    final double isError = recalInfo.getErrorFraction(eventType, offset);

                    incrementDatumOrPutIfNecessary(qualityScoreTable, qual, isError, keys[0], keys[1], eventIndex);

                    for (int i = 2; i < covariates.length; i++) {
                        if (keys[i] < 0)
                            continue;

                        incrementDatumOrPutIfNecessary(tables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
                    }
                }
            }
        }
    }

    /**
     * creates a datum object with one observation and one or zero error
     *
     * @param reportedQual  the quality score reported by the instrument for this base
     * @param isError       whether or not the observation is an error
     * @return a new RecalDatum object with the observation and the error
     */
    protected RecalDatum createDatumObject(final byte reportedQual, final double isError) {
        return new RecalDatum(1, isError, reportedQual);
    }

    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to updateDataForRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public void finalizeData() {
        if ( finalized ) throw new IllegalStateException("FinalizeData() has already been called");

        // merge all of the thread-local tables
        finalRecalibrationTables = mergeThreadLocalRecalibrationTables();

        final NestedIntegerArray<RecalDatum> byReadGroupTable = finalRecalibrationTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> byQualTable = finalRecalibrationTables.getQualityScoreTable();

        // iterate over all values in the qual table
        for ( NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[2];
            final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
            final RecalDatum qualDatum = leaf.value;

            if ( rgDatum == null ) {
                // create a copy of qualDatum, and initialize byReadGroup table with it
                byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            } else {
                // combine the qual datum with the existing datum in the byReadGroup table
                rgDatum.combine(qualDatum);
            }
        }

        finalized = true;
    }

    /**
     * Merge all of the thread local recalibration tables into a single one.
     *
     * Reuses one of the recalibration tables to hold the merged table, so this function can only be
     * called once in the engine.
     *
     * @return the merged recalibration table
     */
    @Requires("! finalized")
    private RecalibrationTables mergeThreadLocalRecalibrationTables() {
        if ( recalibrationTablesList.isEmpty() ) throw new IllegalStateException("recalibration tables list is empty");

        RecalibrationTables merged = null;
        for ( final RecalibrationTables table : recalibrationTablesList ) {
            if ( merged == null )
                // fast path -- if there's only only one table, so just make it the merged one
                merged = table;
            else {
                merged.combine(table);
            }
        }

        return merged;
    }

    /**
     * Get the final recalibration tables, after finalizeData() has been called
     *
     * This returns the finalized recalibration table collected by this engine.
     *
     * It is an error to call this function before finalizeData has been called
     *
     * @return the finalized recalibration table collected by this engine
     */
    public RecalibrationTables getFinalRecalibrationTables() {
        if ( ! finalized ) throw new IllegalStateException("Cannot get final recalibration tables until finalizeData() has been called");
        return finalRecalibrationTables;
    }

    /**
     * Increments the RecalDatum at the specified position in the specified table, or put a new item there
     * if there isn't already one.
     *
     * Does this in a thread-safe way WITHOUT being synchronized: relies on the behavior of NestedIntegerArray.put()
     * to return false if another thread inserts a new item at our position in the middle of our put operation.
     *
     * @param table the table that holds/will hold our item
     * @param qual qual for this event
     * @param isError error value for this event
     * @param keys location in table of our item
     */
    protected void incrementDatumOrPutIfNecessary( final NestedIntegerArray<RecalDatum> table,
                                                   final byte qual,
                                                   final double isError,
                                                   final int... keys ) {
        final RecalDatum existingDatum = table.get(keys);

        if ( existingDatum == null ) {
            // No existing item, try to put a new one
            if ( ! table.put(createDatumObject(qual, isError), keys) ) {
                // Failed to put a new item because another thread came along and put an item here first.
                // Get the newly-put item and increment it (item is guaranteed to exist at this point)
                table.get(keys).increment(1L, isError);
            }
        }
        else {
            // Easy case: already an item here, so increment it
            existingDatum.increment(1L, isError);
        }
    }
}

package org.broadinstitute.sting.gatk.walkers.bqsr;

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

import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.classloader.PublicPackageSource;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class StandardRecalibrationEngine implements RecalibrationEngine, PublicPackageSource {

    protected Covariate[] covariates;
    protected RecalibrationTables recalibrationTables;

    public void initialize(final Covariate[] covariates, final RecalibrationTables recalibrationTables) {
        this.covariates = covariates.clone();
        this.recalibrationTables = recalibrationTables;
    }

    @Override
    public void updateDataForRead( final GATKSAMRecord read, final boolean[] skip, final double[] snpErrors, final double[] insertionErrors, final double[] deletionErrors ) {
        for( int offset = 0; offset < read.getReadBases().length; offset++ ) {
            if( !skip[offset] ) {
                final ReadCovariates readCovariates = covariateKeySetFrom(read);

                final byte qual = read.getBaseQualities()[offset];
                final double isError = snpErrors[offset];

                final int[] keys = readCovariates.getKeySet(offset, EventType.BASE_SUBSTITUTION);
                final int eventIndex = EventType.BASE_SUBSTITUTION.index;

                combineDatumOrPutIfNecessary(recalibrationTables.getReadGroupTable(), qual, isError, keys[0], eventIndex);

                incrementDatumOrPutIfNecessary(recalibrationTables.getQualityScoreTable(), qual, isError, keys[0], keys[1], eventIndex);

                for (int i = 2; i < covariates.length; i++) {
                    if (keys[i] < 0)
                        continue;

                    incrementDatumOrPutIfNecessary(recalibrationTables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
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
     * Get the covariate key set from a read
     *
     * @param read the read
     * @return the covariate keysets for this read
     */
    protected ReadCovariates covariateKeySetFrom(GATKSAMRecord read) {
        return (ReadCovariates) read.getTemporaryAttribute(BaseRecalibrator.COVARS_ATTRIBUTE);
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
    protected void incrementDatumOrPutIfNecessary( final NestedIntegerArray<RecalDatum> table, final byte qual, final double isError, final int... keys ) {
        final RecalDatum existingDatum = table.get(keys);

        if ( existingDatum == null ) {
            // No existing item, try to put a new one
            if ( ! table.put(createDatumObject(qual, isError), keys) ) {
                // Failed to put a new item because another thread came along and put an item here first.
                // Get the newly-put item and increment it (item is guaranteed to exist at this point)
                table.get(keys).increment(1.0, isError);
            }
        }
        else {
            // Easy case: already an item here, so increment it
            existingDatum.increment(1.0, isError);
        }
    }

    /**
     * Combines the RecalDatum at the specified position in the specified table with a new RecalDatum, or put a
     * new item there if there isn't already one.
     *
     * Does this in a thread-safe way WITHOUT being synchronized: relies on the behavior of NestedIntegerArray.put()
     * to return false if another thread inserts a new item at our position in the middle of our put operation.
     *
     * @param table the table that holds/will hold our item
     * @param qual qual for this event
     * @param isError error value for this event
     * @param keys location in table of our item
     */
    protected void combineDatumOrPutIfNecessary( final NestedIntegerArray<RecalDatum> table, final byte qual, final double isError, final int... keys ) {
        final RecalDatum existingDatum = table.get(keys);
        final RecalDatum newDatum = createDatumObject(qual, isError);

        if ( existingDatum == null ) {
            // No existing item, try to put a new one
            if ( ! table.put(newDatum, keys) ) {
                // Failed to put a new item because another thread came along and put an item here first.
                // Get the newly-put item and combine it with our item (item is guaranteed to exist at this point)
                table.get(keys).combine(newDatum);
            }
        }
        else {
            // Easy case: already an item here, so combine it with our item
            existingDatum.combine(newDatum);
        }
    }
}

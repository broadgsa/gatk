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

package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.Covariate;
import org.broadinstitute.sting.gatk.walkers.bqsr.EventType;
import org.broadinstitute.sting.gatk.walkers.bqsr.RecalDatum;
import org.broadinstitute.sting.utils.collections.IntegerIndexedNestedHashMap;

/**
 * Utility class to facilitate on-the-fly base quality score recalibration.
 *
 * User: ebanks
 * Date: 6/20/12
 */

public class RecalibrationTables {

    public enum TableType {
        READ_GROUP_TABLE(0),
        QUALITY_SCORE_TABLE(1),
        OPTIONAL_COVARIATE_TABLES_START(2);

        public final int index;

        private TableType(final int index) {
            this.index = index;
        }
    }

    private final IntegerIndexedNestedHashMap[] tables;

    public RecalibrationTables(final Covariate[] covariates) {
        this(covariates, covariates[TableType.READ_GROUP_TABLE.index].maximumKeyValue() + 1);
    }

    public RecalibrationTables(final Covariate[] covariates, final int numReadGroups) {
        tables = new IntegerIndexedNestedHashMap[covariates.length];

        final int qualDimension = covariates[TableType.QUALITY_SCORE_TABLE.index].maximumKeyValue() + 1;
        final int eventDimension = EventType.values().length;

        tables[TableType.READ_GROUP_TABLE.index] = new IntegerIndexedNestedHashMap<RecalDatum>(numReadGroups, eventDimension);
        tables[TableType.QUALITY_SCORE_TABLE.index] = new IntegerIndexedNestedHashMap<RecalDatum>(numReadGroups, qualDimension, eventDimension);
        for (int i = TableType.OPTIONAL_COVARIATE_TABLES_START.index; i < covariates.length; i++)
            tables[i] = new IntegerIndexedNestedHashMap<RecalDatum>(numReadGroups, qualDimension, covariates[i].maximumKeyValue()+1, eventDimension);
    }

    public IntegerIndexedNestedHashMap<RecalDatum> getTable(final TableType type) {
        return (IntegerIndexedNestedHashMap<RecalDatum>)tables[type.index];
    }

    public IntegerIndexedNestedHashMap<RecalDatum> getTable(final int index) {
        return (IntegerIndexedNestedHashMap<RecalDatum>)tables[index];
    }

    public int numTables() {
        return tables.length;
    }
}

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

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.recalibration.covariates.*;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class RecalibrationTablesUnitTest extends BaseTest {
    @Test
    public void basicTest() {
        final Covariate[] covariates = RecalibrationTestUtils.makeInitializedStandardCovariates();
        final int numReadGroups = 6;
        final RecalibrationTables tables = new RecalibrationTables(covariates, numReadGroups);

        final Covariate qualCov = covariates[1];
        final Covariate cycleCov = covariates[2];
        final Covariate contextCov = covariates[3];

        Assert.assertEquals(tables.numTables(), covariates.length);

        Assert.assertNotNull(tables.getReadGroupTable());
        Assert.assertEquals(tables.getReadGroupTable(), tables.getTable(RecalibrationTables.TableType.READ_GROUP_TABLE.index));
        testDimensions(tables.getReadGroupTable(), numReadGroups);

        Assert.assertNotNull(tables.getQualityScoreTable());
        Assert.assertEquals(tables.getQualityScoreTable(), tables.getTable(RecalibrationTables.TableType.QUALITY_SCORE_TABLE.index));
        testDimensions(tables.getQualityScoreTable(), numReadGroups, qualCov.maximumKeyValue() + 1);

        Assert.assertNotNull(tables.getTable(2));
        testDimensions(tables.getTable(2), numReadGroups, qualCov.maximumKeyValue() + 1, cycleCov.maximumKeyValue() + 1);

        Assert.assertNotNull(tables.getTable(3));
        testDimensions(tables.getTable(3), numReadGroups, qualCov.maximumKeyValue() + 1, contextCov.maximumKeyValue() + 1);
    }

    private void testDimensions(final NestedIntegerArray<RecalDatum> table, final int ... dimensions) {
        final int[] dim = new int[dimensions.length+1];
        System.arraycopy(dimensions, 0, dim, 0, dimensions.length);
        dim[dimensions.length] = EventType.values().length;
        Assert.assertEquals(table.getDimensions().length, dim.length);

        for ( int i = 0; i < dim.length; i++ ) {
            Assert.assertEquals(table.getDimensions()[i], dim[i], "Table dimensions not expected at dim " + i);
        }
    }

    @Test
    public void basicMakeQualityScoreTable() {
        final Covariate[] covariates = RecalibrationTestUtils.makeInitializedStandardCovariates();
        final int numReadGroups = 6;
        final RecalibrationTables tables = new RecalibrationTables(covariates, numReadGroups);

        final Covariate qualCov = covariates[1];
        final NestedIntegerArray<RecalDatum> copy = tables.makeQualityScoreTable();
        testDimensions(copy, numReadGroups, qualCov.maximumKeyValue()+1);
        Assert.assertEquals(copy.getAllValues().size(), 0);
    }
}

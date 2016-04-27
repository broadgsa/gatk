/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import org.broadinstitute.gatk.engine.recalibration.covariates.Covariate;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.collections.NestedIntegerArray;
import org.broadinstitute.gatk.utils.recalibration.EventType;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class RecalibrationTablesUnitTest extends BaseTest {
    private RecalibrationTables tables;
    private Covariate[] covariates;
    private int numReadGroups = 6;
    final byte qualByte = 1;
    final List<Integer> combineStates = Arrays.asList(0, 1, 2);

    @BeforeMethod
    private void makeTables() {
        covariates = RecalibrationTestUtils.makeInitializedStandardCovariates();
        tables = new RecalibrationTables(covariates, numReadGroups);
        fillTable(tables);
    }

    private void fillTable(final RecalibrationTables tables) {
        for ( int iterations = 0; iterations < 10; iterations++ ) {
            for ( final EventType et : EventType.values() ) {
                for ( final int rg : combineStates) {
                    final double error = rg % 2 == 0 ? 1 : 0;
                    RecalUtils.incrementDatumOrPutIfNecessary(tables.getReadGroupTable(), qualByte, error, rg, et.ordinal());
                    for ( final int qual : combineStates) {
                        RecalUtils.incrementDatumOrPutIfNecessary(tables.getQualityScoreTable(), qualByte, error, rg, qual, et.ordinal());
                        for ( final int cycle : combineStates)
                            RecalUtils.incrementDatumOrPutIfNecessary(tables.getTable(2), qualByte, error, rg, qual, cycle, et.ordinal());
                        for ( final int context : combineStates)
                            RecalUtils.incrementDatumOrPutIfNecessary(tables.getTable(3), qualByte, error, rg, qual, context, et.ordinal());
                    }
                }
            }
        }
    }

    @Test
    public void basicTest() {
        final Covariate qualCov = covariates[1];
        final Covariate cycleCov = covariates[2];
        final Covariate contextCov = covariates[3];

        Assert.assertEquals(tables.numTables(), covariates.length);

        Assert.assertNotNull(tables.getReadGroupTable());
        Assert.assertEquals(tables.getReadGroupTable(), tables.getTable(RecalibrationTables.TableType.READ_GROUP_TABLE.ordinal()));
        testDimensions(tables.getReadGroupTable(), numReadGroups);

        Assert.assertNotNull(tables.getQualityScoreTable());
        Assert.assertEquals(tables.getQualityScoreTable(), tables.getTable(RecalibrationTables.TableType.QUALITY_SCORE_TABLE.ordinal()));
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
        final Covariate qualCov = covariates[1];
        final NestedIntegerArray<RecalDatum> copy = tables.makeQualityScoreTable();
        testDimensions(copy, numReadGroups, qualCov.maximumKeyValue()+1);
        Assert.assertEquals(copy.getAllValues().size(), 0);
    }

    @Test
    public void testCombine1() {
        final RecalibrationTables merged = new RecalibrationTables(covariates, numReadGroups);
        fillTable(merged);

        merged.combine(tables);

        for ( int i = 0; i < tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations() * 2);
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches() * 2);
            }
        }
    }

    @Test
    public void testCombineEmptyOther() {
        final RecalibrationTables merged = new RecalibrationTables(covariates, numReadGroups);

        merged.combine(tables);

        for ( int i = 0; i < tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations());
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches());
            }
        }
    }

    @Test
    public void testCombinePartial() {
        final RecalibrationTables merged = new RecalibrationTables(covariates, numReadGroups);
        for ( final int rg : combineStates) {
            RecalUtils.incrementDatumOrPutIfNecessary(merged.getTable(3), qualByte, 1, rg, 0, 0, 0);
        }

        merged.combine(tables);
        for ( int i = 0; i < tables.numTables(); i++ ) {
            NestedIntegerArray<RecalDatum> table = tables.getTable(i);
            NestedIntegerArray<RecalDatum> mergedTable = merged.getTable(i);

            Assert.assertEquals(table.getAllLeaves().size(), mergedTable.getAllLeaves().size());
            for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves() ) {
                final RecalDatum mergedValue = mergedTable.get(leaf.keys);
                Assert.assertNotNull(mergedValue);

                final int delta = i == 3 && leaf.keys[1] == 0 && leaf.keys[2] == 0 && leaf.keys[3] == 0 ? 1 : 0;
                Assert.assertEquals(mergedValue.getNumObservations(), leaf.value.getNumObservations() + delta);
                Assert.assertEquals(mergedValue.getNumMismatches(), leaf.value.getNumMismatches() + delta);
            }
        }
    }
}

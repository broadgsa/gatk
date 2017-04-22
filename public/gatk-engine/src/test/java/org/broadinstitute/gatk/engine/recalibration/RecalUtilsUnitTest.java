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

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.NestedIntegerArray;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public final class RecalUtilsUnitTest extends BaseTest {
    private class Row {
        int rg, qual, ne, no;

        private Row(final Row copy) {
            this(copy.rg, copy.qual, copy.ne, copy.no);
        }

        private Row(int rg, int qual, int ne, int no) {
            this.rg = rg;
            this.qual = qual;
            this.ne = ne;
            this.no = no;
        }

        @Override
        public String toString() {
            return "Row{" +
                    "" + rg +
                    ", " + qual +
                    ", " + ne +
                    ", " + no +
                    '}';
        }
    }

    @DataProvider(name = "CombineTablesProvider")
    public Object[][] createCombineTablesProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Row> rows = new ArrayList<Row>();
        for ( final int rg : Arrays.asList(0, 1) ) {
            for ( final int qual : Arrays.asList(0, 1) ) {
                rows.add(new Row(rg, qual, 1, 10));
            }
        }

        logger.warn("Number of rows " + rows.size());

        List<List<Row>> permutations = new LinkedList<List<Row>>();
        permutations.addAll(Utils.makePermutations(rows, 1, false));
        permutations.addAll(Utils.makePermutations(rows, 2, false));
        permutations.addAll(Utils.makePermutations(rows, 3, false));

        // adding 1 row to 2
        for ( final List<Row> table1 : permutations ) {
            for ( final Row table2 : rows ) {
                tests.add(new Object[]{table1, Arrays.asList(table2)});
            }
        }

        // adding 2 rows to 1
        for ( final List<Row> table1 : permutations ) {
            for ( final Row table2 : rows ) {
                tests.add(new Object[]{Arrays.asList(table2), table1});
            }
        }

        for ( final List<Row> table1 : permutations ) {
            for ( final List<Row> table2 : permutations ) {
                tests.add(new Object[]{table1, table2});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CombineTablesProvider")
    public void testCombineTables(final List<Row> table1, final List<Row> table2) {
        final NestedIntegerArray<RecalDatum> nia1 = makeTable(table1);
        final NestedIntegerArray<RecalDatum> nia2 = makeTable(table2);
        final List<Row> expectedRows = makeExpected(table1, table2);
        final NestedIntegerArray<RecalDatum> expected = makeTable(expectedRows);
        RecalUtils.combineTables(nia1, nia2);

        Assert.assertEquals(nia1.getDimensions(), expected.getDimensions());
        Assert.assertEquals(nia1.getAllValues().size(), expected.getAllValues().size());

        for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : expected.getAllLeaves() ) {
            final RecalDatum actual = nia1.get(leaf.keys);
            Assert.assertEquals(actual.getNumMismatches(), leaf.value.getNumMismatches());
            Assert.assertEquals(actual.getNumObservations(), leaf.value.getNumObservations());
        }
    }

    public List<Row> makeExpected(final List<Row> table1, final List<Row> table2) {
        final List<Row> combined = new LinkedList<Row>();
        for ( final Row t1 : table1 ) combined.add(new Row(t1));
        for ( final Row t2 : table2 ) {
            combine(combined, t2);
        }
        return combined;
    }

    private void combine(final List<Row> combined, final Row row) {
        for ( final Row c : combined ) {
            if ( c.rg == row.rg && c.qual == row.qual ) {
                c.ne += row.ne;
                c.no += row.no;
                return;
            }
        }

        combined.add(new Row(row));
    }

    public NestedIntegerArray<RecalDatum> makeTable(final List<Row> rows) {
        final NestedIntegerArray<RecalDatum> x = new NestedIntegerArray<RecalDatum>(3, 3);
        for ( final Row r : rows )
            x.put(new RecalDatum((long)r.no, (double)r.ne, (byte)10), r.rg, r.qual);
        return x;
    }
}

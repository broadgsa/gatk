/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.report;

import org.broadinstitute.sting.utils.collections.Pair;

import java.util.*;

/**
 * Tracks a linked list of GATKReportColumn in order by name.
 */
public class GATKReportColumns extends LinkedHashMap<String, GATKReportColumn> implements Iterable<GATKReportColumn> {
    private List<String> columnNames = new ArrayList<String>();

    /**
     * Returns the column by index
     * @param i the index
     * @return The column
     */
    public GATKReportColumn getByIndex(int i) {
        return get(columnNames.get(i));
    }

    @Override
    public GATKReportColumn remove(Object key) {
        columnNames.remove(key);
        return super.remove(key);
    }

    @Override
    public GATKReportColumn put(String key, GATKReportColumn value) {
        columnNames.add(key);
        return super.put(key, value);
    }

    @Override
    public Iterator<GATKReportColumn> iterator() {
        return new Iterator<GATKReportColumn>() {
            int offset = 0;
            public boolean hasNext() { return offset < columnNames.size() ; }
            public GATKReportColumn next() { return getByIndex(offset++); }
            public void remove() { throw new UnsupportedOperationException("Cannot remove from a GATKReportColumn iterator"); }
        };
    }
}

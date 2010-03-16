/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.utils.report;

import java.util.*;


/**
 * 
 * @author aaron 
 * 
 * Class TableContainer
 *
 * represents a table of data, and formats according to column count.  Do not use,
 * needs improvement.  This is a stand-in for something better.
 */
public class TableContainer {
    private List<Object> tableEntries;
    private int columns = 2;  // default to two

    public TableContainer(Collection<Object> tbl, int columns) {
        tableEntries = new ArrayList<Object>();
        if (tbl instanceof List || tbl instanceof Vector || tbl instanceof Set)
            tableEntries.addAll(tbl);
        if (tbl instanceof Map)
            for (Object key : ((Map)tbl).keySet()) {
                tableEntries.add(key);
                tableEntries.add(((Map) tbl).get(key));
            }
        this.columns = columns;         
    }

    public List<List<Object>> toRows() {
        ArrayList<List<Object>> list = new ArrayList<List<Object>>();
        List<Object> currentRow = new ArrayList<Object>();
        for (Object obj : tableEntries) {
            if (currentRow.size() >= columns) {
                list.add(currentRow);
                currentRow = new ArrayList<Object>();
            }
            currentRow.add(obj);
        }
        list.add(currentRow);
        return list;
    }

}

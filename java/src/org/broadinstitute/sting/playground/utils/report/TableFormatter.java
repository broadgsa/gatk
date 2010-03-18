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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class TableFormatter
 *
 * a helper class for formatting a table of data.  Given a collection, it formats it into rows and columns.
 * If it's a Map, values are emited as key, value, key, value, so it may make sense to have colums set to 2.
 *
 * TODO: Needs some cleanup
 *
 */
public class TableFormatter {
    public Object tableEntries;                 // data storage
    private int columns = 2;                    // default to two columns of data, overridden in the constructor


    public TableFormatter(Collection<Object> tbl, int columns) {
        tableEntries = tbl;
        this.columns = columns;
    }

    public TableFormatter(Object tbl) {
        tableEntries = tbl;        
    }

    /**
     * convert the underlying object to a list of lists, for convenient output
     * @return a list of list of objects (a list for each row, then a list of row members)
     */
    public List<List<Object>> toRows() {
        if (tableEntries instanceof Collection)
            return collectionToTable(((Collection)tableEntries));
        else if (tableEntries instanceof Map)
            return mapToTable(((Map)tableEntries));
        else
            throw new UnsupportedOperationException("underlying collection must be an instance of type Collection or of type Map");
    }

    private List<List<Object>> mapToTable(Map collection) {
        ArrayList<List<Object>> list = new ArrayList<List<Object>>();        
        for (Object obj : collection.keySet()) {
            List<Object> currentRow = new ArrayList<Object>();
            currentRow.add(obj);
            Object data = collection.get(obj);
            currentRow.add(data.toString());

            list.add(currentRow);
        }

        return list;
    }

    private List<List<Object>> collectionToTable(Collection collection) {
        ArrayList<List<Object>> list = new ArrayList<List<Object>>();
        List<Object> currentRow = new ArrayList<Object>();
        for (Object obj : collection) {
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

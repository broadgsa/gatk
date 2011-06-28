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

package org.broadinstitute.sting.gatk.refdata.features.annotator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.broad.tribble.Feature;

/**
 * This class represents a single record in an AnnotatorInputTable.
 */
public class AnnotatorInputTableFeature implements Feature {

    private ArrayList<String> columnNames;
    private HashMap<String, String> columnValues; //maps colum names to column values

    private String chr;
    private int start;
    private int end;
    private String strRep = null;

    /**
     * Constructor.
     * @param chr    The chromosome name.
     * @param start  The start position
     * @param end    The end position
     */
    public AnnotatorInputTableFeature(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }


    /**
     * Constructor.
     * @param columnNames  The column names as parsed out of the file header.
     */
    public AnnotatorInputTableFeature(ArrayList<String> columnNames) {
        this.columnNames = columnNames;
        this.columnValues = new HashMap<String, String>();
    }



    /**
     * @return the list of column names from the file header.
     */
    public ArrayList<String> getHeader() {
        return columnNames;
    }


    /**
     * Returns the value of the given column.
     *
     * @param columnName The column name as it appears in the file header.
     * @return The value
     */
    public String getColumnValue(final String columnName) {
        return columnValues.get(columnName);
    }


    public boolean containsColumnName(final String columnName) {
        return columnValues.containsKey(columnName);
    }


    /**
     * Sets the value for the given column.
     *
     * @param columnName The column name as it appears in the file header.
     * @param value The value
     * @return The existing value associated with the columnName, if there is one.
     */
    protected String putColumnValue(final String columnName, final String value) {
        return columnValues.put(columnName, value);
    }

    /**
     * @return all values in this line, hashed by their column names.
     */
    public Map<String,String> getColumnValues() {
        return Collections.unmodifiableMap(columnValues);
    }


    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    protected void setChr(String chr) {
        this.chr = chr;
    }

    protected void setStart(int start) {
        this.start = start;
    }

    protected void setEnd(int end) {
        this.end = end;
    }

    @Override
    public String toString() {
        if ( strRep == null ) {
            StringBuilder sb = new StringBuilder();

            for(String columnName : columnNames ) {
                if ( sb.length() == 0 )
                    sb.append("[");
                else
                    sb.append(", ");
                sb.append(columnName + "=" + columnValues.get(columnName));
            }
            sb.append("]");

            strRep = sb.toString();
        }

        return strRep;
    }
}

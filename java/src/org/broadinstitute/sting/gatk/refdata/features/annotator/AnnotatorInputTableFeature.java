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
import java.util.Set;
import java.util.Map.Entry;

import org.broad.tribble.Feature;

/**
 * This class represents a single record in an AnnotatorInputTable.
 */
public class AnnotatorInputTableFeature implements Feature {

    private ArrayList<String> columnNames;
    private HashMap<String, String> columnValues;

    private String chr;
    private int start;
    private int end;


    // Temporary attributes were added to make it easier to implement certain
    // optimizations for RODs that span an interval. For example, if a Walker
    // needs to do a time-consuming computation on data from a ROD, it would normally
    // have to repeat this computation every time its map(..) method is called.
    // If a ROD spans an interval, the Walker's map(..) method will be called for every position in ROD.
    // However, many computations (including validation and parsing) are done per ROD rather than
    // per position. Therefore, substantial optimizations are possible if the result
    // of the first computation is cached and reused on subsequent map(..) calls.
    // Temporary attributes provide a convenient place to store these results,
    // freeing the Walkers from having to maintain their own ROD -> result hashmaps.
    private Map<Object, Object> temporaryAttributes;




    /**
     * Constructor.
     * @param columnNames  The column names as parsed out of the file header.
     */
    public AnnotatorInputTableFeature(ArrayList<String> columnNames) {
        this.columnNames = columnNames;
        this.columnValues = new HashMap<String, String>();
    }



    /**
     * Returns the list of column names from the file header.
     * @return
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
    public String getColumnValue(final Object columnName) {
        return columnValues.get(columnName);
    }


    public boolean containsColumnName(final Object columnName) {
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
     * Returns all values in this line, hashed by their column names.
     *
     * @return
     */
    public Map<String,String> getColumnValues() {
        return Collections.unmodifiableMap(columnValues);
    }


    /**
     * Returns the entry set of all column name-value pairs.
     *
     * @return
     */
    public Set<Entry<String, String>> getEntrySet() {

        return columnValues.entrySet();
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


    /**
     * Checks whether an attribute has been set for the given key.
     *
     * Temporary attributes make it easier to implement certain
     * optimizations for RODs that span an interval. For example, if a Walker
     * needs to do a time-consuming computation on data from a ROD, it would normally
     * have to repeat this computation every time its map(..) method is called.
     * If a ROD spans an interval, the Walker's map(..) method will be called for every position in ROD.
     * However, many computations (including validation and parsing) are done per ROD rather than
     * per position. Therefore, substantial optimizations are possible if the result
     * of the first computation is cached and reused on subsequent map(..) calls.
     * Temporary attributes provide a convenient place to store these results,
     * freeing the Walkers from having to maintain their own ROD -> result hashmaps.
     *
     * @param key key
     * @return True if an attribute has been set for this key.
     */
    public boolean containsTemporaryAttribute(Object key) {
        if(temporaryAttributes != null) {
            return temporaryAttributes.containsKey(key);
        }
        return false;
    }

    /**
     * Sets the key to the given value, replacing any previous value. The previous
     * value is returned.
     *
     * Temporary attributes make it easier to implement certain
     * optimizations for RODs that span an interval. For example, if a Walker
     * needs to do a time-consuming computation on data from a ROD, it would normally
     * have to repeat this computation every time its map(..) method is called.
     * If a ROD spans an interval, the Walker's map(..) method will be called for every position in ROD.
     * However, many computations (including validation and parsing) are done per ROD rather than
     * per position. Therefore, substantial optimizations are possible if the result
     * of the first computation is cached and reused on subsequent map(..) calls.
     * Temporary attributes provide a convenient place to store these results,
     * freeing the Walkers from having to maintain their own ROD -> result hashmaps.
     *
     * @param key    key
     * @param value  value
     * @return attribute
     */
    public Object setTemporaryAttribute(Object key, Object value) {
        if(temporaryAttributes == null) {
            temporaryAttributes = new HashMap<Object, Object>();
        }
        return temporaryAttributes.put(key, value);
    }

    /**
     * Looks up the value associated with the given key.
     *
     * Temporary attributes make it easier to implement certain
     * optimizations for RODs that span an interval. For example, if a Walker
     * needs to do a time-consuming computation on data from a ROD, it would normally
     * have to repeat this computation every time its map(..) method is called.
     * If a ROD spans an interval, the Walker's map(..) method will be called for every position in ROD.
     * However, many computations (including validation and parsing) are done per ROD rather than
     * per position. Therefore, substantial optimizations are possible if the result
     * of the first computation is cached and reused on subsequent map(..) calls.
     * Temporary attributes provide a convenient place to store these results,
     * freeing the Walkers from having to maintain their own ROD -> result hashmaps.
     *
     * @param key key
     * @return The value, or null.
     */
    public Object getTemporaryAttribute(Object key) {
        if(temporaryAttributes != null) {
            return temporaryAttributes.get(key);
        }
        return null;
    }

    /**
     * Removes the attribute that has the given key.
     *
     * Temporary attributes make it easier to implement certain
     * optimizations for RODs that span an interval. For example, if a Walker
     * needs to do a time-consuming computation on data from a ROD, it would normally
     * have to repeat this computation every time its map(..) method is called.
     * If a ROD spans an interval, the Walker's map(..) method will be called for every position in ROD.
     * However, many computations (including validation and parsing) are done per ROD rather than
     * per position. Therefore, substantial optimizations are possible if the result
     * of the first computation is cached and reused on subsequent map(..) calls.
     * Temporary attributes provide a convenient place to store these results,
     * freeing the Walkers from having to maintain their own ROD -> result hashmaps.
     *
     * @param key key
     * @return The value that was associated with this key, or null.
     */
    public Object removeTemporaryAttribute(Object key) {
         if(temporaryAttributes != null) {
             return temporaryAttributes.remove(key);
         }
         return null;
    }





}

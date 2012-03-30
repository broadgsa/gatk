/*
 * Copyright (c) 2012, The Broad Institute
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

import org.apache.commons.lang.ObjectUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GATKReportTable {
    /**
     * REGEX that matches any table with an invalid name
     */
    public static final String INVALID_TABLE_NAME_REGEX = "[^a-zA-Z0-9_\\-\\.]";
    private static final String GATKTABLE_HEADER_PREFIX = "#:GATKTable";
    private static final String SEPARATOR = ":";
    private static final String ENDLINE = ":;";

    private String tableName;
    private String tableDescription;


    private String primaryKeyName;
    private Collection<Object> primaryKeyColumn;
    private boolean primaryKeyDisplay;
    private boolean sortByPrimaryKey = true;

    private GATKReportColumns columns;

    private static final String COULD_NOT_READ_HEADER = "Could not read the header of this file -- ";
    private static final String COULD_NOT_READ_COLUMN_NAMES = "Could not read the column names of this file -- ";
    private static final String COULD_NOT_READ_DATA_LINE = "Could not read a data line of this table -- ";
    private static final String COULD_NOT_READ_EMPTY_LINE = "Could not read the last empty line of this table -- ";
    private static final String OLD_GATK_TABLE_VERSION = "We no longer support older versions of the GATK Tables";
    
    private static final String NUMBER_CONVERSION_EXCEPTION = "String is a number but is not a long or a double: ";

    public GATKReportTable(BufferedReader reader, GATKReportVersion version) {
        int counter = 0;

        switch (version) {
            case V1_0:
                int nHeaders = 2;
                String[] tableHeaders = new String[nHeaders];
    
                // Read in the headers
                for (int i = 0; i < nHeaders; i++) {
                    try {
                        tableHeaders[i] = reader.readLine();
                    } catch (IOException e) {
                        throw new ReviewedStingException(COULD_NOT_READ_HEADER + e.getMessage());
                    }
                }
                String[] tableData = tableHeaders[0].split(":");
                String[] userData = tableHeaders[1].split(":");
    
                // Fill in the fields
                tableName = userData[2];
                tableDescription = (userData.length <= 3) ? "" : userData[3];                                           // table may have no description! (and that's okay)
                primaryKeyDisplay = Boolean.parseBoolean(tableData[2]);
                columns = new GATKReportColumns();
    
                int nColumns = Integer.parseInt(tableData[3]);
                int nRows = Integer.parseInt(tableData[4]);
    
    
                // Read column names
                String columnLine;
                try {
                columnLine = reader.readLine();
                } catch (IOException e) {
                    throw new ReviewedStingException(COULD_NOT_READ_COLUMN_NAMES);
                }
    
                List<Integer> columnStarts = TextFormattingUtils.getWordStarts(columnLine);
                String[] columnNames = TextFormattingUtils.splitFixedWidth(columnLine, columnStarts);
    
                if (primaryKeyDisplay) {
                    addPrimaryKey(columnNames[0]);
    
                } else {
                    sortByPrimaryKey = true;
                    addPrimaryKey("id", false);
                    counter = 1;
                }
                // Put in columns using the format string from the header
                for (int i = 0; i < nColumns; i++) {
                    String format = tableData[5 + i];
                    if (primaryKeyDisplay)
                        addColumn(columnNames[i + 1], true, format);
                    else
                        addColumn(columnNames[i], true, format);
                }
    
                for (int i = 0; i < nRows; i++) {
                    // read line
                    String dataLine;
                    try {
                        dataLine = reader.readLine(); 
                    } catch (IOException e) {
                        throw new ReviewedStingException(COULD_NOT_READ_DATA_LINE + e.getMessage());
                    }
                    List<String> lineSplits = Arrays.asList(TextFormattingUtils.splitFixedWidth(dataLine, columnStarts));
    
                    for (int columnIndex = 0; columnIndex < nColumns; columnIndex++) {
    
                        //Input all the remaining values
                        GATKReportDataType type = getColumns().getByIndex(columnIndex).getDataType();
    
                        if (primaryKeyDisplay) {
                            String columnName = columnNames[columnIndex + 1];
                            String primaryKey = lineSplits.get(0);
                            set(primaryKey, columnName, type.Parse(lineSplits.get(columnIndex + 1)));
                        } else {
                            String columnName = columnNames[columnIndex];
                            set(counter, columnName, type.Parse(lineSplits.get(columnIndex)));
                        }
    
                    }
                    counter++;
                }
    
    
                try {
                    reader.readLine();
                } catch (IOException e) {
                    throw new ReviewedStingException(COULD_NOT_READ_EMPTY_LINE + e.getMessage());
                }  
            break;
            
            default: 
                throw new ReviewedStingException(OLD_GATK_TABLE_VERSION);
        }
    }


    /**
     * Verifies that a table or column name has only alphanumeric characters - no spaces or special characters allowed
     *
     * @param name the name of the table or column
     * @return true if the name is valid, false if otherwise
     */
    private boolean isValidName(String name) {
        Pattern p = Pattern.compile(INVALID_TABLE_NAME_REGEX);
        Matcher m = p.matcher(name);

        return !m.find();
    }

    /**
     * Verifies that a table or column name has only alphanumeric characters - no spaces or special characters allowed
     *
     * @param description the name of the table or column
     * @return true if the name is valid, false if otherwise
     */
    private boolean isValidDescription(String description) {
        Pattern p = Pattern.compile("\\r|\\n");
        Matcher m = p.matcher(description);

        return !m.find();
    }

    /**
     * Construct a new GATK report table with the specified name and description
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     */
    public GATKReportTable(String tableName, String tableDescription) {
        this(tableName, tableDescription, true);
    }

    /**
     * Construct a new GATK report table with the specified name and description and whether to sort rows by the primary
     * key
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param sortByPrimaryKey whether to sort rows by the primary key (instead of order added)
     */
    public GATKReportTable(String tableName, String tableDescription, boolean sortByPrimaryKey) {
        if (!isValidName(tableName)) {
            throw new ReviewedStingException("Attempted to set a GATKReportTable name of '" + tableName + "'.  GATKReportTable names must be purely alphanumeric - no spaces or special characters are allowed.");
        }

        if (!isValidDescription(tableDescription)) {
            throw new ReviewedStingException("Attempted to set a GATKReportTable description of '" + tableDescription + "'.  GATKReportTable descriptions must not contain newlines.");
        }

        this.tableName = tableName;
        this.tableDescription = tableDescription;
        this.sortByPrimaryKey = sortByPrimaryKey;

        columns = new GATKReportColumns();
    }

    /**
     * Add a primary key column.  This becomes the unique identifier for every column in the table.
     *
     * @param primaryKeyName the name of the primary key column
     */
    public void addPrimaryKey(String primaryKeyName) {
        addPrimaryKey(primaryKeyName, true);
    }

    /**
     * Add an optionally visible primary key column.  This becomes the unique identifier for every column in the table,
     * and will always be printed as the first column.
     *
     * @param primaryKeyName the name of the primary key column
     * @param display        should this primary key be displayed?
     */
    public void addPrimaryKey(String primaryKeyName, boolean display) {
        if (!isValidName(primaryKeyName)) {
            throw new ReviewedStingException("Attempted to set a GATKReportTable primary key name of '" + primaryKeyName + "'.  GATKReportTable primary key names must be purely alphanumeric - no spaces or special characters are allowed.");
        }

        this.primaryKeyName = primaryKeyName;

        primaryKeyColumn = sortByPrimaryKey ? new TreeSet<Object>() : new LinkedList<Object>();
        primaryKeyDisplay = display;
    }

    /**
     * Returns the first primary key matching the dotted column values.
     * Ex: dbsnp.eval.called.all.novel.all
     *
     * @param dottedColumnValues Period concatenated values.
     * @return The first primary key matching the column values or throws an exception.
     */
    public Object getPrimaryKeyByData(String dottedColumnValues) {
        Object key = findPrimaryKey(dottedColumnValues);
        if (key == null)
            throw new ReviewedStingException("Attempted to get non-existent GATKReportTable key for values: " + dottedColumnValues);
        return key;
    }

    /**
     * Returns true if there is at least on row with the dotted column values.
     * Ex: dbsnp.eval.called.all.novel.all
     *
     * @param dottedColumnValues Period concatenated values.
     * @return true if there is at least one row matching the columns.
     */
    public boolean containsPrimaryKey(String dottedColumnValues) {
        return findPrimaryKey(dottedColumnValues) != null;
    }

    /**
     * Returns the first primary key matching the dotted column values.
     * Ex: dbsnp.eval.called.all.novel.all
     *
     * @param dottedColumnValues Period concatenated values.
     * @return The first primary key matching the column values or null.
     */
    private Object findPrimaryKey(String dottedColumnValues) {
        return findPrimaryKey(dottedColumnValues.split("\\."));
    }

    /**
     * Returns the first primary key matching the column values.
     * Ex: new String[] { "dbsnp", "eval", "called", "all", "novel", "all" }
     *
     * @param columnValues column values.
     * @return The first primary key matching the column values.
     */
    private Object findPrimaryKey(Object[] columnValues) {
        for (Object primaryKey : primaryKeyColumn) {
            boolean matching = true;
            for (int i = 0; matching && i < columnValues.length; i++) {
                matching = ObjectUtils.equals(columnValues[i], get(primaryKey, i + 1));
            }
            if (matching)
                return primaryKey;
        }
        return null;
    }

    /**
     * Add a column to the report and specify the default value that should be supplied if a given position in the table
     * is never explicitly set.
     *
     * @param columnName   the name of the column
     * @param defaultValue the default value for the column
     */
    public void addColumn(String columnName, Object defaultValue) {
        addColumn(columnName, defaultValue, true);
    }

    /**
     * Add a column to the report, specify the default column value, and specify whether the column should be displayed
     * in the final output (useful when intermediate columns are necessary for later calculations, but are not required
     * to be in the output file.
     *
     * @param columnName   the name of the column
     * @param defaultValue the default value of the column
     * @param display      if true - the column will be displayed; if false - the column will be hidden
     */
    public void addColumn(String columnName, Object defaultValue, boolean display) {
        addColumn(columnName, defaultValue, display, "");
    }

    /**
     * Add a column to the report, specify the default column value, and specify whether the column should be displayed
     * in the final output (useful when intermediate columns are necessary for later calculations, but are not required
     * to be in the output file.
     *
     * @param columnName   the name of the column
     * @param defaultValue the default value of the column
     * @param format       the format string used to display data
     */
    public void addColumn(String columnName, Object defaultValue, String format) {
        addColumn(columnName, defaultValue, true, format);
    }

    /**
     * Add a column to the report, specify whether the column should be displayed in the final output (useful when
     * intermediate columns are necessary for later calculations, but are not required to be in the output file), and the
     * format string used to display the data.
     *
     * @param columnName the name of the column
     * @param display    if true - the column will be displayed; if false - the column will be hidden
     * @param format     the format string used to display data
     */
    public void addColumn(String columnName, boolean display, String format) {
        addColumn(columnName, null, display, format);
    }

    /**
     * Add a column to the report, specify the default column value, whether the column should be displayed in the final
     * output (useful when intermediate columns are necessary for later calculations, but are not required to be in the
     * output file), and the format string used to display the data.
     *
     * @param columnName   the name of the column
     * @param defaultValue the default value of a blank cell
     * @param display      if true - the column will be displayed; if false - the column will be hidden
     * @param format       the format string used to display data
     */
    public void addColumn(String columnName, Object defaultValue, boolean display, String format) {
        if (!isValidName(columnName)) {
            throw new ReviewedStingException("Attempted to set a GATKReportTable column name of '" + columnName + "'.  GATKReportTable column names must be purely alphanumeric - no spaces or special characters are allowed.");
        }
        columns.put(columnName, new GATKReportColumn(columnName, defaultValue, display, format));
    }

    /**
     * Check if the requested element exists, and if not, create it.
     *
     * @param primaryKey the primary key value
     * @param columnName the name of the column
     */
    private void verifyEntry(Object primaryKey, String columnName) {
        if (!columns.containsKey(columnName)) {
            throw new ReviewedStingException("Attempted to access column '" + columnName + "' that does not exist in table '" + tableName + "'.");
        }

        primaryKeyColumn.add(primaryKey);

        if (!columns.get(columnName).containsKey(primaryKey)) {
            columns.get(columnName).initialize(primaryKey);
        }
    }

    public boolean containsKey(Object primaryKey) {
        return primaryKeyColumn.contains(primaryKey);
    }

    public Collection<Object> getPrimaryKeys() {
        return Collections.unmodifiableCollection(primaryKeyColumn);
    }

    /**
     * Set the value for a given position in the table
     *
     * @param primaryKey the primary key value
     * @param columnName the name of the column
     * @param value      the value to set
     */
    public void set(Object primaryKey, String columnName, Object value) {
        verifyEntry(primaryKey, columnName);
        GATKReportColumn column = columns.get(columnName);
        //todo -- Check if value is of same type as column

        // We do not accept internal null values
        if (value == null)
            value = "null";

        // This code below is bs. Why am do I have to conform to bad code
        // Below is some code to convert a string into its appropriate type.
        
        // I second Roger's rant!

        // If we got a string but the column is not a String type
        Object newValue = null;
        if (value instanceof String && !column.getDataType().equals(GATKReportDataType.String)) {
            // Integer case
            if (column.getDataType().equals(GATKReportDataType.Integer)) {
                try {
                    newValue = Long.parseLong((String) value);
                } catch (Exception e) {
                    /** do nothing */
                }
            }
            if (column.getDataType().equals(GATKReportDataType.Decimal)) {
                try {
                    newValue = Double.parseDouble((String) value);
                } catch (Exception e) {
                    /** do nothing */
                }
            }
            if (column.getDataType().equals(GATKReportDataType.Character) && ((String) value).length() == 1) {
                newValue = ((String) value).charAt(0);
            }
        }

        if (newValue != null)
            value = newValue;

        // todo -- Types have to be more flexible. For example, %d should accept Integers, Shorts and Bytes.
        if (column.getDataType().equals(GATKReportDataType.fromObject(value)) || column.getDataType().equals(GATKReportDataType.Unknown) )
            columns.get(columnName).put(primaryKey, value);
        else
            throw new ReviewedStingException(String.format("Tried to add an object of type: %s to a column of type: %s", GATKReportDataType.fromObject(value).name(), column.getDataType().name()));
    }

    /**
     * Get a value from the given position in the table
     *
     * @param primaryKey the primary key value
     * @param columnName the name of the column
     * @return the value stored at the specified position in the table
     */
    public Object get(Object primaryKey, String columnName) {
        verifyEntry(primaryKey, columnName);

        return columns.get(columnName).get(primaryKey);
    }

    /**
     * Get a value from the given position in the table
     *
     * @param primaryKey  the primary key value
     * @param columnIndex the index of the column
     * @return the value stored at the specified position in the table
     */
    private Object get(Object primaryKey, int columnIndex) {
        return columns.getByIndex(columnIndex).get(primaryKey);
    }

    /**
     * Increment an element in the table.  This implementation is awful - a functor would probably be better.
     *
     * @param primaryKey the primary key value
     * @param columnName the name of the column
     */
    public void increment(Object primaryKey, String columnName) {
        Object oldValue = get(primaryKey, columnName);
        Object newValue;

        if (oldValue instanceof Byte) {
            newValue = ((Byte) oldValue) + 1;
        } else if (oldValue instanceof Short) {
            newValue = ((Short) oldValue) + 1;
        } else if (oldValue instanceof Integer) {
            newValue = ((Integer) oldValue) + 1;
        } else if (oldValue instanceof Long) {
            newValue = ((Long) oldValue) + 1L;
        } else if (oldValue instanceof Float) {
            newValue = ((Float) oldValue) + 1.0f;
        } else if (oldValue instanceof Double) {
            newValue = ((Double) oldValue) + 1.0d;
        } else {
            throw new UnsupportedOperationException("Can't increment value: object " + oldValue + " is not an instance of one of the numerical Java primitive wrapper classes (Byte, Short, Integer, Long, Float, Double)");
        }

        set(primaryKey, columnName, newValue);
    }

    /**
     * Decrement an element in the table.  This implementation is awful - a functor would probably be better.
     *
     * @param primaryKey the primary key value
     * @param columnName the name of the column
     */
    public void decrement(Object primaryKey, String columnName) {
        Object oldValue = get(primaryKey, columnName);
        Object newValue;

        if (oldValue instanceof Byte) {
            newValue = ((Byte) oldValue) - 1;
        } else if (oldValue instanceof Short) {
            newValue = ((Short) oldValue) - 1;
        } else if (oldValue instanceof Integer) {
            newValue = ((Integer) oldValue) - 1;
        } else if (oldValue instanceof Long) {
            newValue = ((Long) oldValue) - 1L;
        } else if (oldValue instanceof Float) {
            newValue = ((Float) oldValue) - 1.0f;
        } else if (oldValue instanceof Double) {
            newValue = ((Double) oldValue) - 1.0d;
        } else {
            throw new UnsupportedOperationException("Can't decrement value: object " + oldValue + " is not an instance of one of the numerical Java primitive wrapper classes (Byte, Short, Integer, Long, Float, Double)");
        }

        set(primaryKey, columnName, newValue);
    }

    /**
     * Add the specified value to an element in the table
     *
     * @param primaryKey the primary key value
     * @param columnName the name of the column
     * @param valueToAdd the value to add
     */
    public void add(Object primaryKey, String columnName, Object valueToAdd) {
        Object oldValue = get(primaryKey, columnName);
        Object newValue;

        if (oldValue instanceof Byte) {
            newValue = ((Byte) oldValue) + ((Byte) valueToAdd);
        } else if (oldValue instanceof Short) {
            newValue = ((Short) oldValue) + ((Short) valueToAdd);
        } else if (oldValue instanceof Integer) {
            newValue = ((Integer) oldValue) + ((Integer) valueToAdd);
        } else if (oldValue instanceof Long) {
            newValue = ((Long) oldValue) + ((Long) valueToAdd);
        } else if (oldValue instanceof Float) {
            newValue = ((Float) oldValue) + ((Float) valueToAdd);
        } else if (oldValue instanceof Double) {
            newValue = ((Double) oldValue) + ((Double) valueToAdd);
        } else {
            throw new UnsupportedOperationException("Can't add values: object " + oldValue + " is not an instance of one of the numerical Java primitive wrapper classes (Byte, Short, Integer, Long, Float, Double)");
        }

        set(primaryKey, columnName, newValue);
    }

    /**
     * Subtract the specified value from an element in the table
     *
     * @param primaryKey      the primary key value
     * @param columnName      the name of the column
     * @param valueToSubtract the value to subtract
     */
    public void subtract(Object primaryKey, String columnName, Object valueToSubtract) {
        Object oldValue = get(primaryKey, columnName);
        Object newValue;

        if (oldValue instanceof Byte) {
            newValue = ((Byte) oldValue) - ((Byte) valueToSubtract);
        } else if (oldValue instanceof Short) {
            newValue = ((Short) oldValue) - ((Short) valueToSubtract);
        } else if (oldValue instanceof Integer) {
            newValue = ((Integer) oldValue) - ((Integer) valueToSubtract);
        } else if (oldValue instanceof Long) {
            newValue = ((Long) oldValue) - ((Long) valueToSubtract);
        } else if (oldValue instanceof Float) {
            newValue = ((Float) oldValue) - ((Float) valueToSubtract);
        } else if (oldValue instanceof Double) {
            newValue = ((Double) oldValue) - ((Double) valueToSubtract);
        } else {
            throw new UnsupportedOperationException("Can't subtract values: object " + oldValue + " is not an instance of one of the numerical Java primitive wrapper classes (Byte, Short, Integer, Long, Float, Double)");
        }

        set(primaryKey, columnName, newValue);
    }

    /**
     * Multiply the specified value to an element in the table
     *
     * @param primaryKey      the primary key value
     * @param columnName      the name of the column
     * @param valueToMultiply the value to multiply by
     */
    public void multiply(Object primaryKey, String columnName, Object valueToMultiply) {
        Object oldValue = get(primaryKey, columnName);
        Object newValue;

        if (oldValue instanceof Byte) {
            newValue = ((Byte) oldValue) * ((Byte) valueToMultiply);
        } else if (oldValue instanceof Short) {
            newValue = ((Short) oldValue) * ((Short) valueToMultiply);
        } else if (oldValue instanceof Integer) {
            newValue = ((Integer) oldValue) * ((Integer) valueToMultiply);
        } else if (oldValue instanceof Long) {
            newValue = ((Long) oldValue) * ((Long) valueToMultiply);
        } else if (oldValue instanceof Float) {
            newValue = ((Float) oldValue) * ((Float) valueToMultiply);
        } else if (oldValue instanceof Double) {
            newValue = ((Double) oldValue) * ((Double) valueToMultiply);
        } else {
            throw new UnsupportedOperationException("Can't multiply values: object " + oldValue + " is not an instance of one of the numerical Java primitive wrapper classes (Byte, Short, Integer, Long, Float, Double)");
        }

        set(primaryKey, columnName, newValue);
    }

    /**
     * Divide the specified value from an element in the table
     *
     * @param primaryKey    the primary key value
     * @param columnName    the name of the column
     * @param valueToDivide the value to divide by
     */
    public void divide(Object primaryKey, String columnName, Object valueToDivide) {
        Object oldValue = get(primaryKey, columnName);
        Object newValue;

        if (oldValue instanceof Byte) {
            newValue = ((Byte) oldValue) * ((Byte) valueToDivide);
        } else if (oldValue instanceof Short) {
            newValue = ((Short) oldValue) * ((Short) valueToDivide);
        } else if (oldValue instanceof Integer) {
            newValue = ((Integer) oldValue) * ((Integer) valueToDivide);
        } else if (oldValue instanceof Long) {
            newValue = ((Long) oldValue) * ((Long) valueToDivide);
        } else if (oldValue instanceof Float) {
            newValue = ((Float) oldValue) * ((Float) valueToDivide);
        } else if (oldValue instanceof Double) {
            newValue = ((Double) oldValue) * ((Double) valueToDivide);
        } else {
            throw new UnsupportedOperationException("Can't divide values: object " + oldValue + " is not an instance of one of the numerical Java primitive wrapper classes (Byte, Short, Integer, Long, Float, Double)");
        }

        set(primaryKey, columnName, newValue);
    }

    /**
     * Add two columns to each other and set the results to a third column
     *
     * @param columnToSet the column that should hold the results
     * @param augend      the column that shall be the augend
     * @param addend      the column that shall be the addend
     */
    public void addColumns(String columnToSet, String augend, String addend) {
        for (Object primaryKey : primaryKeyColumn) {
            Number firstColumnValue = (Number) get(primaryKey, augend);
            Number secondColumnValue = (Number) get(primaryKey, addend);

            Double value = firstColumnValue.doubleValue() + secondColumnValue.doubleValue();

            set(primaryKey, columnToSet, value);
        }
    }

    /**
     * Subtract one column from another and set the results to a third column
     *
     * @param columnToSet the column that should hold the results
     * @param minuend     the column that shall be the minuend (the a in a - b)
     * @param subtrahend  the column that shall be the subtrahend (the b in a - b)
     */
    public void subtractColumns(String columnToSet, String minuend, String subtrahend) {
        for (Object primaryKey : primaryKeyColumn) {
            Number firstColumnValue = (Number) get(primaryKey, minuend);
            Number secondColumnValue = (Number) get(primaryKey, subtrahend);

            Double value = firstColumnValue.doubleValue() - secondColumnValue.doubleValue();

            set(primaryKey, columnToSet, value);
        }
    }

    /**
     * Multiply two columns by each other and set the results to a third column
     *
     * @param columnToSet  the column that should hold the results
     * @param multiplier   the column that shall be the multiplier
     * @param multiplicand the column that shall be the multiplicand
     */
    public void multiplyColumns(String columnToSet, String multiplier, String multiplicand) {
        for (Object primaryKey : primaryKeyColumn) {
            Number firstColumnValue = (Number) get(primaryKey, multiplier);
            Number secondColumnValue = (Number) get(primaryKey, multiplicand);

            Double value = firstColumnValue.doubleValue() * secondColumnValue.doubleValue();

            set(primaryKey, columnToSet, value);
        }
    }

    /**
     * Divide two columns by each other and set the results to a third column
     *
     * @param columnToSet       the column that should hold the results
     * @param numeratorColumn   the column that shall be the numerator
     * @param denominatorColumn the column that shall be the denominator
     */
    public void divideColumns(String columnToSet, String numeratorColumn, String denominatorColumn) {
        for (Object primaryKey : primaryKeyColumn) {
            Number firstColumnValue = (Number) get(primaryKey, numeratorColumn);
            Number secondColumnValue = (Number) get(primaryKey, denominatorColumn);

            Double value = firstColumnValue.doubleValue() / secondColumnValue.doubleValue();

            set(primaryKey, columnToSet, value);
        }
    }

    /**
     * Return the print width of the primary key column
     *
     * @return the width of the primary key column
     */
    int getPrimaryKeyColumnWidth() {
        int maxWidth = getPrimaryKeyName().length();

        for (Object primaryKey : primaryKeyColumn) {
            int width = primaryKey.toString().length();

            if (width > maxWidth) {
                maxWidth = width;
            }
        }

        return maxWidth;
    }

    /**
     * Write the table to the PrintStream, formatted nicely to be human-readable, AWK-able, and R-friendly.
     *
     * @param out the PrintStream to which the table should be written
     */
     void write(PrintStream out) {

        /*
         * Table header:
         * #:GATKTable:nColumns:nRows:(DataType for each column):;
         * #:GATKTable:TableName:Description :;
         * key   colA  colB
         * row1  xxxx  xxxxx
        */

        // Get the column widths for everything
        String primaryKeyFormat = "%-" + getPrimaryKeyColumnWidth() + "s";

        // Emit the table definition
        String formatHeader = String.format(GATKTABLE_HEADER_PREFIX + ":%b:%d:%d", primaryKeyDisplay, getColumns().size(), getNumRows());
        // Add all the formats for all the columns
        for (GATKReportColumn column : getColumns()) {
            if (column.isDisplayable())
                formatHeader += (SEPARATOR + column.getFormat());
        }
        out.println(formatHeader + ENDLINE);
        out.printf(GATKTABLE_HEADER_PREFIX + ":%s:%s\n", tableName, tableDescription);

        //out.printf("#:GATKTable:%s:%s", Algorithm);


        // Emit the table header, taking into account the padding requirement if the primary key is a hidden column
        boolean needsPadding = false;
        if (primaryKeyDisplay) {
            out.printf(primaryKeyFormat, getPrimaryKeyName());
            needsPadding = true;
        }

        for (String columnName : columns.keySet()) {
            if (columns.get(columnName).isDisplayable()) {
                if (needsPadding) {
                    out.printf("  ");
                }
                out.printf(columns.get(columnName).getColumnFormat().getNameFormat(), columnName);

                needsPadding = true;
            }
        }

        out.printf("%n");

        // Emit the table body
        for (final Object primaryKey : primaryKeyColumn) {
            needsPadding = false;
            if (primaryKeyDisplay) {
                out.printf(primaryKeyFormat, primaryKey);
                needsPadding = true;
            }

            for (final Map.Entry<String, GATKReportColumn> entry : columns.entrySet()) {
                final GATKReportColumn column = entry.getValue();
                if (column.isDisplayable()) {
                    if (needsPadding) {
                        out.print("  ");
                    }

                    final String value = column.getStringValue(primaryKey);
                    out.printf(column.getColumnFormat().getValueFormat(), value);

                    needsPadding = true;
                }
            }

            out.println();
        }

         out.println();
    }

    public int getNumRows() {
        return primaryKeyColumn.size();
    }

    public String getTableName() {
        return tableName;
    }

    public String getTableDescription() {
        return tableDescription;
    }

    public GATKReportColumns getColumns() {
        return columns;
    }

    /**
     * Combines two compatible GATK report tables. This is the general function which will call the different algorithms
     * necessary to gather the tables. Every column's combine algorithm is read and treated accordingly.
     *
     * @param input Another GATK table
     */
    void combineWith(GATKReportTable input) {
        /*
         * This function is different from addRowsFrom because we will add the ability to sum,average, etc rows
         * TODO: Add other combining algorithms
         */

        // Make sure the columns match AND the Primary Key
        if (input.getColumns().keySet().equals(this.getColumns().keySet()) &&
                input.getPrimaryKeyName().equals(this.getPrimaryKeyName())) {
            this.addRowsFrom(input);
        } else
            throw new ReviewedStingException("Failed to combine GATKReportTable, columns don't match!");
    }

    /**
     * A gather algorithm that simply takes the rows from the argument, and adds them to the current table. This is the
     * default gather algorithm.
     *
     * @param input Another GATK table to add rows from.
     */
    private void addRowsFrom(GATKReportTable input) {
        // add column by column

        // For every column
        for (String columnKey : input.getColumns().keySet()) {
            GATKReportColumn current = this.getColumns().get(columnKey);
            GATKReportColumn toAdd = input.getColumns().get(columnKey);
            // We want to take the current column and add all the values from input

            // The column is a map of values <Key, Value>
            for (Object rowKey : toAdd.keySet()) {
                // We add every value from toAdd to the current
                if (!current.containsKey(rowKey)) {
                    this.set(rowKey, columnKey, toAdd.get(rowKey));
                    //System.out.printf("Putting row with PK: %s \n", rowKey);
                } else {
                    this.set(rowKey, columnKey, toAdd.get(rowKey));

                    System.out.printf("OVERWRITING Row with PK: %s \n", rowKey);
                }
            }
        }

    }

    public String getPrimaryKeyName() {
        return primaryKeyName;
    }

    /**
     * Returns whether or not the two tables have the same format including columns and everything in between. This does
     * not check if the data inside is the same. This is the check to see if the two tables are gatherable or
     * reduceable
     *
     * @param table another GATK table
     * @return true if the the tables are gatherable
     */
    public boolean isSameFormat(GATKReportTable table) {
        //Should we add the sortByPrimaryKey as a check?

        return columns.isSameFormat(table.columns) &&
               (primaryKeyDisplay == table.primaryKeyDisplay && primaryKeyName.equals(table.primaryKeyName) &&
               tableName.equals(table.tableName) &&
               tableDescription.equals(table.tableDescription));
    }

    /**
     * Checks that the tables are exactly the same.
     *
     * @param table another GATK report
     * @return true if all field in the reports, tables, and columns are equal.
     */
    public boolean equals(GATKReportTable table) {
        return isSameFormat(table) &&
               (columns.equals(table.columns) &&
               primaryKeyColumn.equals(table.primaryKeyColumn) &&
               sortByPrimaryKey == table.sortByPrimaryKey);

    }
}

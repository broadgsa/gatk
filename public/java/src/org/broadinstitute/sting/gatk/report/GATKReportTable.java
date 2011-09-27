package org.broadinstitute.sting.gatk.report;

import org.apache.commons.lang.ObjectUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A data structure that allows data to be collected over the course of a walker's computation, then have that data
 * written to a PrintStream such that it's human-readable, AWK-able, and R-friendly (given that you load it using the
 * GATKReport loader module).
 *
 * The goal of this object is to use the same data structure for both accumulating data during a walker's computation
 * and emitting that data to a file for easy analysis in R (or any other program/language that can take in a table of
 * results).  Thus, all of the infrastructure below is designed simply to make printing the following as easy as
 * possible:
 *
 * ##:GATKReport.v0.1 ErrorRatePerCycle : The error rate per sequenced position in the reads
 * cycle  errorrate.61PA8.7         qualavg.61PA8.7
 * 0      0.007451835696110506      25.474613284804366
 * 1      0.002362777171937477      29.844949954504095
 * 2      9.087604507451836E-4      32.87590975254731
 * 3      5.452562704471102E-4      34.498999090081895
 * 4      9.087604507451836E-4      35.14831665150137
 * 5      5.452562704471102E-4      36.07223435225619
 * 6      5.452562704471102E-4      36.1217248908297
 * 7      5.452562704471102E-4      36.1910480349345
 * 8      5.452562704471102E-4      36.00345705967977
 *
 * Here, we have a GATKReport table - a well-formatted, easy to read representation of some tabular data.  Every single
 * table has this same GATKReport.v0.1 header, which permits multiple files from different sources to be cat-ed
 * together, which makes it very easy to pull tables from different programs into R via a single file.
 *
 * ------------
 * Definitions:
 *
 * Table info:
 *   The first line, structured as
 *     ##:<report version> <table name> : <table description>
 *
 * Table header:
 *   The second line, specifying a unique name for each column in the table.
 *
 *   The first column mentioned in the table header is the "primary key" column - a column that provides the unique
 *   identifier for each row in the table.  Once this column is created, any element in the table can be referenced by
 *   the row-column coordinate, i.e. "primary key"-"column name" coordinate.
 *
 *   When a column is added to a table, a default value must be specified (usually 0).  This is the initial value for
 *   an element in a column.  This permits operations like increment() and decrement() to work properly on columns that
 *   are effectively counters for a particular event.
 *
 *   Finally, the display property for each column can be set during column creation.  This is useful when a given
 *   column stores an intermediate result that will be used later on, perhaps to calculate the value of another column.
 *   In these cases, it's obviously necessary to store the value required for further computation, but it's not
 *   necessary to actually print the intermediate column.
 *
 * Table body:
 *   The values of the table itself.
 *
 * ---------------
 * Implementation:
 *
 * The implementation of this table has two components:
 *   1. A TreeSet<Object> that stores all the values ever specified for the primary key.  Any get() operation that
 *      refers to an element where the primary key object does not exist will result in its implicit creation.  I
 *      haven't yet decided if this is a good idea...
 *
 *   2. A HashMap<String, GATKReportColumn> that stores a mapping from column name to column contents.  Each
 *      GATKReportColumn is effectively a map (in fact, GATKReportColumn extends TreeMap<Object, Object>) between
 *      primary key and the column value.  This means that, given N columns, the primary key information is stored
 *      N+1 times.  This is obviously wasteful and can likely be handled much more elegantly in future implementations.
 *
 * ------------------------------
 * Element and column operations:
 *
 * In addition to simply getting and setting values, this object also permits some simple operations to be applied to
 * individual elements or to whole columns.  For instance, an element can be easily incremented without the hassle of
 * calling get(), incrementing the obtained value by 1, and then calling set() with the new value.  Also, some vector
 * operations are supported.  For instance, two whole columns can be divided and have the result be set to a third
 * column.  This is especially useful when aggregating counts in two intermediate columns that will eventually need to
 * be manipulated row-by-row to compute the final column.
 *
 * Note: I've made no attempt whatsoever to make these operations efficient.  Right now, some of the methods check the
 * type of the stored object using an instanceof call and attempt to do the right thing.  Others cast the contents of
 * the cell to a Number, call the Number.toDouble() method and compute a result.  This is clearly not the ideal design,
 * but at least the prototype contained herein works.
 *
 * @author Kiran Garimella
 * @author Khalid Shakir
 */
public class GATKReportTable {
    /** REGEX that matches any table with an invalid name */
    public final static String INVALID_TABLE_NAME_REGEX = "[^a-zA-Z0-9_\\-\\.]";
    private static final GATKReportVersion LATEST_REPORT_VERSION = GATKReportVersion.V0_2;
    private String tableName;
    private String tableDescription;
    private GATKReportVersion version = LATEST_REPORT_VERSION;

    private String primaryKeyName;
    private Collection<Object> primaryKeyColumn;
    private boolean primaryKeyDisplay;
    private boolean sortByPrimaryKey = true;

    private GATKReportColumns columns;

    /**
     * Verifies that a table or column name has only alphanumeric characters - no spaces or special characters allowed
     *
     * @param name  the name of the table or column
     * @return  true if the name is valid, false if otherwise
     */
    private boolean isValidName(String name) {
        Pattern p = Pattern.compile(INVALID_TABLE_NAME_REGEX);
        Matcher m = p.matcher(name);

        return !m.find();
    }

    /**
     * Verifies that a table or column name has only alphanumeric characters - no spaces or special characters allowed
     *
     * @param description  the name of the table or column
     * @return  true if the name is valid, false if otherwise
     */
    private boolean isValidDescription(String description) {
        Pattern p = Pattern.compile("\\r|\\n");
        Matcher m = p.matcher(description);

        return !m.find();
    }

    /**
     * Construct a new GATK report table with the specified name and description
     *
     * @param tableName  the name of the table
     * @param tableDescription  the description of the table
     */
    public GATKReportTable(String tableName, String tableDescription) {
        this(tableName, tableDescription, true);
    }

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

    public GATKReportVersion getVersion() {
        return version;
    }

    protected void setVersion(GATKReportVersion version) {
        this.version = version;
    }

    /**
     * Add a primary key column.  This becomes the unique identifier for every column in the table.
     *
     * @param primaryKeyName  the name of the primary key column
     */
    public void addPrimaryKey(String primaryKeyName) {
        addPrimaryKey(primaryKeyName, true);
    }

    /**
     * Add an optionally visible primary key column.  This becomes the unique identifier for every column in the table, and will always be printed as the first column.
     *
     * @param primaryKeyName  the name of the primary key column
     * @param display should this primary key be displayed?
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
     * @param dottedColumnValues Period concatenated values.
     * @return The first primary key matching the column values or throws an exception.
     */
    public Object getPrimaryKey(String dottedColumnValues) {
        Object key = findPrimaryKey(dottedColumnValues);
        if (key == null)
            throw new ReviewedStingException("Attempted to get non-existent GATKReportTable key for values: " + dottedColumnValues);
        return key;
    }

    /**
     * Returns true if there is at least on row with the dotted column values.
     * Ex: dbsnp.eval.called.all.novel.all
     * @param dottedColumnValues Period concatenated values.
     * @return true if there is at least one row matching the columns.
     */
    public boolean containsPrimaryKey(String dottedColumnValues) {
        return findPrimaryKey(dottedColumnValues) != null;
    }

    /**
     * Returns the first primary key matching the dotted column values.
     * Ex: dbsnp.eval.called.all.novel.all
     * @param dottedColumnValues Period concatenated values.
     * @return The first primary key matching the column values or null.
     */
    private Object findPrimaryKey(String dottedColumnValues) {
        return findPrimaryKey(dottedColumnValues.split("\\."));
    }

    /**
     * Returns the first primary key matching the column values.
     * Ex: new String[] { "dbsnp", "eval", "called", "all", "novel", "all" }
     * @param columnValues column values.
     * @return The first primary key matching the column values.
     */
    private Object findPrimaryKey(Object[] columnValues) {
        for (Object primaryKey : primaryKeyColumn) {
            boolean matching = true;
            for (int i = 0; matching && i < columnValues.length; i++) {
                matching = ObjectUtils.equals(columnValues[i], get(primaryKey, i+1));
            }
            if (matching)
                return primaryKey;
        }
        return null;
    }

    /**
     * Add a column to the report and specify the default value that should be supplied if a given position in the table is never explicitly set.
     *
     * @param columnName  the name of the column
     * @param defaultValue  the default value for the column
     */
    public void addColumn(String columnName, Object defaultValue) {
        if (!isValidName(columnName)) {
            throw new ReviewedStingException("Attempted to set a GATKReportTable column name of '" + columnName + "'.  GATKReportTable column names must be purely alphanumeric - no spaces or special characters are allowed.");
        }

        addColumn(columnName, defaultValue, true);
    }

    /**
     * Add a column to the report, specify the default column value, and specify whether the column should be displayed in the final output (useful when intermediate columns are necessary for later calculations, but are not required to be in the output file.
     *
     * @param columnName  the name of the column
     * @param defaultValue  the default value of the column
     * @param display  if true - the column will be displayed; if false - the column will be hidden
     */
    public void addColumn(String columnName, Object defaultValue, boolean display) {
        columns.put(columnName, new GATKReportColumn(columnName, defaultValue, display));
    }

    /**
     * Check if the requested element exists, and if not, create it.
     *
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
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

    /**
     * Set the value for a given position in the table
     *
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
     * @param value  the value to set
     */
    public void set(Object primaryKey, String columnName, Object value) {
        verifyEntry(primaryKey, columnName);

        columns.get(columnName).put(primaryKey, value);
    }

    /**
     * Get a value from the given position in the table
     *
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
     * @return  the value stored at the specified position in the table
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
     * @return  the value stored at the specified position in the table
     */
    private Object get(Object primaryKey, int columnIndex) {
        return columns.getByIndex(columnIndex).get(primaryKey);
    }

    /**
     * Increment an element in the table.  This implementation is awful - a functor would probably be better.
     *
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
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
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
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
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
     * @param valueToAdd  the value to add
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
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
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
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
     * @param valueToMultiply  the value to multiply by
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
     * @param primaryKey  the primary key value
     * @param columnName  the name of the column
     * @param valueToDivide  the value to divide by
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
     * @param columnToSet  the column that should hold the results
     * @param augend  the column that shall be the augend
     * @param addend  the column that shall be the addend
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
     * @param columnToSet  the column that should hold the results
     * @param minuend  the column that shall be the minuend (the a in a - b)
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
     * @param multiplier  the column that shall be the multiplier
     * @param multiplicand  the column that shall be the multiplicand
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
     * @param columnToSet  the column that should hold the results
     * @param numeratorColumn  the column that shall be the numerator
     * @param denominatorColumn  the column that shall be the denominator
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
     * @return  the width of the primary key column
     */
    public int getPrimaryKeyColumnWidth() {
        int maxWidth = primaryKeyName.length();

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
     * @param out  the PrintStream to which the table should be written
     */
    public void write(PrintStream out) {
        // Get the column widths for everything
        HashMap<String, String> columnWidths = new HashMap<String, String>();
        for (String columnName : columns.keySet()) {
            int width = columns.get(columnName).getColumnWidth();
            String format = "%-" + String.valueOf(width) + "s";

            columnWidths.put(columnName, format);
        }
        String primaryKeyFormat = "%-" + getPrimaryKeyColumnWidth() + "s";

        // Emit the table definition
        out.printf("##:GATKReport.%s %s : %s%n", LATEST_REPORT_VERSION.versionString, tableName, tableDescription);

        // Emit the table header, taking into account the padding requirement if the primary key is a hidden column
        boolean needsPadding = false;
        if (primaryKeyDisplay) {
            out.printf(primaryKeyFormat, primaryKeyName);
            needsPadding = true;
        }

        for (String columnName : columns.keySet()) {
            if (columns.get(columnName).isDisplayable()) {
                if (needsPadding) { out.printf("  "); }
                out.printf(columnWidths.get(columnName), columnName);

                needsPadding = true;
            }
        }

        out.printf("%n");

        // Emit the table body
        for (Object primaryKey : primaryKeyColumn) {
            needsPadding = false;
            if (primaryKeyDisplay) {
                out.printf(primaryKeyFormat, primaryKey);
                needsPadding = true;
            }

            for (String columnName : columns.keySet()) {
                if (columns.get(columnName).isDisplayable()) {
                    if (needsPadding) { out.printf("  "); }
                    String value = columns.get(columnName).getStringValue(primaryKey);
                    out.printf(columnWidths.get(columnName), value);

                    needsPadding = true;
                }
            }

            out.printf("%n");
        }

        // Close the table
        out.printf("%n");
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
}

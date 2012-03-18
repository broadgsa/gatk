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

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Collection;
import java.util.TreeMap;

/**
 * Container class for GATK report tables
 */
public class GATKReport {
    public static final String GATKREPORT_HEADER_PREFIX = "#:GATKReport.";
    public static final GATKReportVersion LATEST_REPORT_VERSION = GATKReportVersion.V1_0;
    private static final String SEPARATOR = ":";
    private GATKReportVersion version = LATEST_REPORT_VERSION;

    private final TreeMap<String, GATKReportTable> tables = new TreeMap<String, GATKReportTable>();

    /**
     * Create a new, empty GATKReport.
     */
    public GATKReport() {
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     *
     * @param filename the path to the file to load
     */
    public GATKReport(String filename) {
        this(new File(filename));
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     *
     * @param file the file to load
     */
    public GATKReport(File file) {
        loadReport(file);
    }

    /**
     * Create a new GATK report from GATK report tables
     * @param tables Any number of tables that you want ot add to the report
     */
    public GATKReport(GATKReportTable... tables) {
        for( GATKReportTable table: tables)
            addTable(table);
    }

    /**
     * Load a GATKReport file from disk
     *
     * @param file the file to load
     */
    private void loadReport(File file) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));

            String reportHeader = reader.readLine();

            // Read the first line for the version and number of tables.
            version = GATKReportVersion.fromHeader(reportHeader);
            if (version.equals(GATKReportVersion.V0_1) ||
                    version.equals(GATKReportVersion.V0_2))
                throw new UserException("The GATK no longer supports reading legacy GATK Reports. Please use v1.0 or newer.");

            int nTables = Integer.parseInt(reportHeader.split(":")[2]);

            // Read each tables according ot the number of tables
            for (int i = 0; i < nTables; i++) {
                addTable(new GATKReportTable(reader, version));

                /*
                if ( !blankLine.equals("") ) {
                    throw new StingException("The GATK Report File is corrupted or not formatted correctly");
                }
                */
            }


        } catch (Exception e) {
            // todo - improve exception handling
            //throw new StingException("Cannot read GATKReport: " + e);
            e.printStackTrace();
        }
    }

    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     */
    public void addTable(String tableName, String tableDescription) {
        addTable(tableName, tableDescription, true);
    }

    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param sortByPrimaryKey whether to sort the rows by the primary key
     */
    public void addTable(String tableName, String tableDescription, boolean sortByPrimaryKey) {
        GATKReportTable table = new GATKReportTable(tableName, tableDescription, sortByPrimaryKey);
        tables.put(tableName, table);
    }

    /**
     * Adds a table, empty or populated, to the report
     *
     * @param table the table to add
     */
    public void addTable(GATKReportTable table) {
        tables.put(table.getTableName(), table);
    }

    /**
     * Return true if table with a given name exists
     *
     * @param tableName the name of the table
     * @return true if the table exists, false otherwise
     */
    public boolean hasTable(String tableName) {
        return tables.containsKey(tableName);
    }

    /**
     * Return a table with a given name
     *
     * @param tableName the name of the table
     * @return the table object
     */
    public GATKReportTable getTable(String tableName) {
        GATKReportTable table = tables.get(tableName);
        if (table == null)
            throw new ReviewedStingException("Table is not in GATKReport: " + tableName);
        return table;
    }

    /**
     * Print all tables contained within this container to a PrintStream
     *
     * @param out the PrintStream to which the tables should be written
     */
    public void print(PrintStream out) {
        out.println(GATKREPORT_HEADER_PREFIX + getVersion().toString() + SEPARATOR + getTables().size());
        for (GATKReportTable table : tables.values()) {
            if (table.getNumRows() > 0) {
                table.write(out);
            }
        }
    }

    public Collection<GATKReportTable> getTables() {
        return tables.values();
    }

    /**
     * This is the main function is charge of gathering the reports. It checks that the reports are compatible and then
     * calls the table atheirng functions.
     *
     * @param input another GATKReport of the same format
     */
    public void combineWith(GATKReport input) {

        if (!this.isSameFormat(input)) {
            throw new ReviewedStingException("Failed to combine GATKReport, format doesn't match!");
        }

        for (String tableName : input.tables.keySet()) {
            tables.get(tableName).combineWith(input.getTable(tableName));
        }

    }

    public GATKReportVersion getVersion() {
        return version;
    }

    /**
     * Returns whether or not the two reports have the same format, from columns, to tables, to reports, and everything
     * in between. This does not check if the data inside is the same. This is the check to see if the two reports are
     * gatherable or reduceable.
     *
     * @param report another GATK report
     * @return true if the the reports are gatherable
     */
    public boolean isSameFormat(GATKReport report) {
        if (!version.equals(report.version)) {
            return false;
        }
        if (!tables.keySet().equals(report.tables.keySet())) {
            return false;
        }
        for (String tableName : tables.keySet()) {
            if (!getTable(tableName).isSameFormat(report.getTable(tableName)))
                return false;
        }
        return true;
    }

    /**
     * Checks that the reports are exactly the same.
     *
     * @param report another GATK report
     * @return true if all field in the reports, tables, and columns are equal.
     */
    public boolean equals(GATKReport report) {
        if (!version.equals(report.version)) {
            return false;
        }
        if (!tables.keySet().equals(report.tables.keySet())) {
            return false;
        }
        for (String tableName : tables.keySet()) {
            if (!getTable(tableName).equals(report.getTable(tableName)))
                return false;
        }
        return true;
    }

    /**
     * The constructor for a simplified GATK Report. Simplified GATK report are designed for reports that do not need
     * the advanced functionality of a full GATK Report.
     * <p/>
     * A simple GATK Report consists of:
     * <p/>
     * - A single table
     * - No primary key ( it is hidden )
     * <p/>
     * Optional:
     * - Only untyped columns. As long as the data is an Object, it will be accepted.
     * - Default column values being empty strings.
     * <p/>
     * Limitations:
     * <p/>
     * - A simple GATK report cannot contain multiple tables.
     * - It cannot contain typed columns, which prevents arithmetic gathering.
     *
     * @param tableName The name of your simple GATK report table
     * @param columns   The names of the columns in your table
     * @return a simplified GATK report
     */
    public static GATKReport newSimpleReport(String tableName, String... columns) {
        GATKReportTable table = new GATKReportTable(tableName, "A simplified GATK table report");
        table.addPrimaryKey("id", false);

        for (String column : columns) {
            table.addColumn(column, "");
        }

        GATKReport output = new GATKReport();
        output.addTable(table);

        return output;
    }

    /**
     * This method provides an efficient way to populate a simplified GATK report. This method will only work on reports
     * that qualify as simplified GATK reports. See the newSimpleReport() constructor for more information.
     *
     * @param values the row of data to be added to the table.
     *               Note: the number of arguments must match the columns in the table.
     */
    public void addRow(Object... values) {
        // Must be a simplified GATK Report
        if (isSimpleReport()) {

            GATKReportTable table = tables.firstEntry().getValue();
            if (table.getColumns().size() != values.length) {
                throw new StingException("The number of arguments in addRow() must match the number of columns in the table");
            }

            int counter = table.getNumRows() + 1;
            int i = 0;

            for (String columnName : table.getColumns().keySet()) {
                table.set(counter, columnName, values[i]);
                i++;
            }

        } else {
            throw new StingException("Cannot add a Row to a non-Simplified GATK Report");
        }


    }

    /**
     * Checks if the GATK report qualifies as a "simple" GATK report
     *
     * @return true is the report is a simplified GATK report
     */
    private boolean isSimpleReport() {
        if (tables.size() != 1)
            return false;

        GATKReportTable table = tables.firstEntry().getValue();

        if (!table.getPrimaryKeyName().equals("id"))
            return false;

        return true;

    }
}

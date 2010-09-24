package org.broadinstitute.sting.gatk.report;

import java.io.PrintStream;
import java.util.TreeMap;

/**
 * Container class for GATK report tables
 */
public class GATKReport {
    private TreeMap<String, GATKReportTable> tables;

    public GATKReport() {
        tables = new TreeMap<String, GATKReportTable>();
    }

    /**
     * Add a new table to the collection
     *
     * @param tableName  the name of the table
     * @param tableDescription  the description of the table
     */
    public void addTable(String tableName, String tableDescription) {
        GATKReportTable table = new GATKReportTable(tableName, tableDescription);
        tables.put(tableName, table);
    }

    /**
     * Return a table with a given name
     *
     * @param tableName  the name of the table
     * @return  the name of the table
     */
    public GATKReportTable getTable(String tableName) {
        return tables.get(tableName);
    }

    /**
     * Print all tables contained within this container to a PrintStream
     *
     * @param out  the PrintStream to which the tables should be written
     */
    public void print(PrintStream out) {
        for (GATKReportTable table : tables.values()) {
            table.write(out);
        }
    }
}

package org.broadinstitute.sting.gatk.report;

import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.*;
import java.util.TreeMap;

/**
 * Container class for GATK report tables
 */
public class GATKReport {
    private TreeMap<String, GATKReportTable> tables;

    /**
     * Create a new, empty GATKReport.
     */
    public GATKReport() {
        tables = new TreeMap<String, GATKReportTable>();
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     * @param filename  the path to the file to load
     */
    public GATKReport(String filename) {
        loadReport(new File(filename));
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     * @param file  the file to load
     */
    public GATKReport(File file) {
        tables = new TreeMap<String, GATKReportTable>();
        loadReport(file);
    }

    /**
     * Load a GATKReport file from disk
     * @param file  the file to load
     */
    private void loadReport(File file) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));

            GATKReportTable table = null;
            String[] header = null;
            int id = 0;

            String line;
            while ( (line = reader.readLine()) != null ) {
                if (line.startsWith("##:GATKReport.v0.1 ")) {
                    line = line.replaceFirst("##:GATKReport.v0.1 ", "");
                    String[] pieces = line.split(" : ");

                    String tableName = pieces[0];
                    String tableDesc = pieces[1];

                    addTable(tableName, tableDesc);
                    table = getTable(tableName);

                    header = null;
                } else if ( line.isEmpty() ) {
                    // do nothing
                } else {
                    if (table != null) {
                        if (header == null) {
                            header = line.split("\\s+");

                            table.addPrimaryKey("id", false);

                            for ( String columnName : header ) {
                                table.addColumn(columnName, "");
                            }

                            id = 0;
                        } else {
                            String[] entries = line.split("\\s+");

                            for (int columnIndex = 0; columnIndex < header.length; columnIndex++) {
                                table.set(id, header[columnIndex], entries[columnIndex]);
                            }

                            id++;
                        }
                    }
                }
            }
        } catch (FileNotFoundException e) {
            throw new StingException("Cannot read GATKReport: " + e);
        } catch (IOException e) {
            throw new StingException("Cannot read GATKReport: " + e);
        }
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
     * Return true if table with a given name exists
     *
     * @param tableName  the name of the table
     * @return true if the table exists, false otherwise
     */
    public boolean hasTable(String tableName) {
        return tables.containsKey(tableName);
    }

    /**
     * Return a table with a given name
     *
     * @param tableName  the name of the table
     * @return  the table object
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
            if (table.getNumRows() > 0) {
                table.write(out);
            }
        }
    }
}

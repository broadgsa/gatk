package org.broadinstitute.sting.gatk.report;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;

import java.io.*;
import java.util.Collection;
import java.util.List;
import java.util.TreeMap;

/**
 * Container class for GATK report tables
 */
public class GATKReport {
    public static final String GATKREPORT_HEADER_PREFIX = "##:GATKReport.v";
    private TreeMap<String, GATKReportTable> tables = new TreeMap<String, GATKReportTable>();

    /**
     * Create a new, empty GATKReport.
     */
    public GATKReport() {
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     * @param filename  the path to the file to load
     */
    public GATKReport(String filename) {
        this(new File(filename));
    }

    /**
     * Create a new GATKReport with the contents of a GATKReport on disk.
     * @param file  the file to load
     */
    public GATKReport(File file) {
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
            GATKReportVersion version = null;
            List<Integer> columnStarts = null;

            String line;
            while ( (line = reader.readLine()) != null ) {

                if (line.startsWith(GATKREPORT_HEADER_PREFIX)) {

                    version = GATKReportVersion.fromHeader(line);

                    line = line.replaceFirst("##:GATKReport." + version.versionString + " ", "");
                    String[] pieces = line.split(" : ");

                    String tableName = pieces[0];
                    String tableDesc = pieces[1];

                    addTable(tableName, tableDesc);
                    table = getTable(tableName);
                    table.setVersion(version);

                    header = null;
                    columnStarts = null;
                } else if ( line.trim().isEmpty() ) {
                    // do nothing
                } else {
                    if (table != null) {

                        String[] splitLine;

                        switch (version) {
                            case V0_1:
                                splitLine = TextFormattingUtils.splitWhiteSpace(line);
                                break;

                            case V0_2:
                                if (header == null) {
                                    columnStarts = TextFormattingUtils.getWordStarts(line);
                                }
                                splitLine = TextFormattingUtils.splitFixedWidth(line, columnStarts);
                                break;

                            default:
                                throw new ReviewedStingException("GATK report version parsing not implemented for: " + line);
                        }

                        if (header == null) {
                            header = splitLine;

                            table.addPrimaryKey("id", false);

                            for ( String columnName : header ) {
                                table.addColumn(columnName, "");
                            }

                            id = 0;
                        } else {
                            for (int columnIndex = 0; columnIndex < header.length; columnIndex++) {
                                table.set(id, header[columnIndex], splitLine[columnIndex]);
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
        addTable(tableName, tableDescription, true);
    }

    public void addTable(String tableName, String tableDescription, boolean sortByPrimaryKey) {
        GATKReportTable table = new GATKReportTable(tableName, tableDescription, sortByPrimaryKey);
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
        GATKReportTable table = tables.get(tableName);
        if (table == null)
            throw new ReviewedStingException("Table is not in GATKReport: " + tableName);
        return table;
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

    public Collection<GATKReportTable> getTables() {
        return tables.values();
    }
}

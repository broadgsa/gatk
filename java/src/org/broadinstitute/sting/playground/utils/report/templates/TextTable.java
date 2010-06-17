package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.utils.StingException;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

/**
 * this class handles generating text table output
 */
class TextTable {
    // make a table of the number of rows
    ArrayList<ArrayList<String>> rows = new ArrayList<ArrayList<String>>();

    // keep a complementary array of the sizes of each column
    ArrayList<Integer> minimumSizes = new ArrayList<Integer>();
    List<String> header = new ArrayList<String>();

    // our default width
    private static final int defaultMinimumSize = 15;

    // our seperator
    private String seperator = "";

    // our default cell entry
    private static final String defaultCell = "";

    private String name;
    private String description;


    /**
     * create a text table, with the:
     * @param name name of the table
     * @param description what the table represents
     * @param header the header fields
     */
    TextTable(String name, String description, List<String> header) {
        this.name = name;
        this.description = description;
        for (int index = 0; index < header.size(); index++)
            determineMinimumColumnWidth(header.get(index).length(),index);
        this.header.addAll(header);
    }

    /**
     * create a text table
     */
    TextTable() {
    }

    /**
     * set a cell of the table
     * @param row, zero based
     * @param column zero based
     * @param entry cell to set
     */
    public void setCell(int row, int column, String entry) {
        if (rows.size() <= row)
            createRow(row,column);
        ArrayList<String> rowData = rows.get(row);
        if (rowData.size() <= column)
            for (int x = rowData.size(); x <= column; x++)
                rowData.add("");
        rows.get(row).set(column,entry);
        determineMinimumColumnWidth(entry.length(),column);
    }

    private void createRow(int row, int column) {
        for (int x = rows.size(); x <= row; x++) {
            ArrayList<String> blank = new ArrayList<String>();
            for (int col = 0; col <= column; col++)
                blank.add(defaultCell);
            rows.add(blank);

        }
    }

    /**
     * write the table to the writer
     * @param writer the writer
     */
    public void toPrintWriter(Writer writer) {
        int index = 0;
        for (String col : header)
            writeToDisk(writer, col, index++, (index == header.size() -1));
        appendEndLine(writer,"\n");
        index = 0;
        for (String col : header)
            writeToDisk(writer, col, index++, (index == header.size() -1));
        appendEndLine(writer,"\n");
        for (ArrayList<String> row : rows) {
            if (row.size() > header.size())
                throw new IllegalStateException("More row-cells in table " + name + " then header columns");
            for (int y = 0; y < header.size(); y++)
                if (row.size() >= y)
                    writeToDisk(writer, "", y, (y == header.size() - 1));
                else
                    writeToDisk(writer, row.get(y), y, (y == header.size() - 1));
        }
        try {
            writer.append("\n");
        } catch (IOException e) {
            throw new StingException("Unable to write to the Writer");
        }
    }

    private void appendEndLine(Writer writer, String str) {
        try {
            writer.append(str);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    /**
     * write to disk a single string
     * @param writer the writer to use
     * @param str the string to write
     * @param y an index of the string, for column width look-up
     * @param lastVal are we the last value
     */
    void writeToDisk(Writer writer, String str, int y, boolean lastVal) {
        try {
            writer.append(String.format("%1$-" + (getMinimumSizes(y) + 3) + "s", str));
            if (y != rows.size() - 1)
                writer.append(seperator);
        } catch (IOException e) {
            throw new StingException("Unable to write to the Writer");
        }
    }

    /**
     * get the minimum size for a column
     * @param column the column index
     * @return the width
     */
    public int getMinimumSizes(int column) {
        if (column >= minimumSizes.size())
            return defaultMinimumSize;
        return minimumSizes.get(column);
    }
    /**
     * determine the minimum column width
     * @param size the size of this string
     * @param position the position (column)
     */
    void determineMinimumColumnWidth(int size, int position) {
        if (minimumSizes.size() <= position)
            for (int x = minimumSizes.size(); x <= position; x++)
                minimumSizes.add(x,0);
        minimumSizes.set(position,(minimumSizes.get(position) > size) ? minimumSizes.get(position) : size);
    }

    public List<String> getHeader() {
        return header;
    }

    public void setHeader(List<String> header) {
        this.header = header;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }
}

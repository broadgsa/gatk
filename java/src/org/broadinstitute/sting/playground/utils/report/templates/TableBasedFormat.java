package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;

import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * an abstract class to share the basics of a table based format; many methods
 * overlap in different output types.
 */
public abstract class TableBasedFormat implements ReportFormat {
    private Map<String, List<Node>> analyses = new HashMap<String, List<Node>>();
    private PrintWriter stream;
    
    @Override
    public void write(Writer writeTo, Node baseNode) {
        getAnalyses(baseNode);
        stream = new PrintWriter(writeTo);
        for (String s : analyses.keySet()) {
            writeAnalysis(analyses.get(s));
            outputTables(analyses.get(s));
        }
    }

    /**
     * break out the analyses by type, given the base node
     * @param baseNode the root node
     */
    private void getAnalyses(Node baseNode) {
        for (Node n : baseNode.getChildren())
            if (!n.tag && n.getComplex()) {
                if (!analyses.containsKey(n.getValue()))
                    analyses.put(n.getValue(),new ArrayList<Node>());
                analyses.get(n.getValue()).add(n);
            }
    }

    /**
     * write the analysis nodes out, only outputting the simple data points (non-table data)
     * @param nodes a list of nodes, of the same analysis type
     */
    private void writeAnalysis(List<Node> nodes) {
        if (nodes.size() < 1 || !nodes.get(0).getName().equals("analysis")) return;
        Node forTitle = nodes.get(0);
        stream.println(niceDivider(80));
        stream.println("Analysis Name:         \t" + forTitle.getValue());
        stream.println("Analysis Description:  \t" + forTitle.getDescription());
        stream.println();

        String header = extractHeaderString(forTitle);
        if (header == null) return; // a null here indicates we don't have any unique columns to display
        stream.println(header);
        stream.println(niceDivider(header.length()));

        for (Node analysis : nodes) {
            String dataString = dataPointNodesToValues(analysis);
            if (dataString.length() > 0 && !dataString.equals("<null>")) {
                stream.print(getTagValues(analysis));
                stream.println(dataString);
            }
        }
        stream.println();
        stream.println();

    }

    /**
     * output the tables: look at list of analysis nodes (all from the same analysis) and output the table
     * @param nodes the list of analysis nodes (of the same underlying type)
     */
    public void outputTables(List<Node> nodes) {
        Map<String,List<String>> tableRows = new HashMap<String,List<String>>();
        Map<String,String> tableHeaders = new HashMap<String,String>();
        for (Node analysis : nodes)
            for (Node n : analysis.getChildren()) {
                if (n.table) {
                    StringBuilder columnBuilder = new StringBuilder();
                    getTagNames(analysis,columnBuilder);
                    for (Node row : n.getChildren()) {
                        StringBuilder rowBuilder = new StringBuilder();
                        rowBuilder.append(getTagValues(analysis));
                        rowBuilder.append(formatColumn(row.getValue()));
                        columnBuilder.append(formatColumn(row.getName()));
                        for (Node column : row.getChildren()) {
                            columnBuilder.append(formatColumn(column.getValue()));
                            if (column.getChildren().size() == 1) {
                                String value =  formatColumn(column.getChildren().iterator().next().getValue());
                                rowBuilder.append(value);
                            }
                        }
                        if (!tableRows.containsKey(n.getValue()))
                            tableRows.put(n.getValue(),new ArrayList<String>());
                        tableRows.get(n.getValue()).add(rowBuilder.toString());
                        if (!tableHeaders.containsKey(n.getValue()))
                            tableHeaders.put(n.getValue(),columnBuilder.toString());
                    }
                }
            }

        // output the tables
        for (String tableName : tableHeaders.keySet()) {
            stream.println("Table Name : " + tableName);
            stream.println();
            stream.println(tableHeaders.get(tableName));
            stream.println(niceDivider(tableHeaders.get(tableName).length()));
            List<String> rows = tableRows.get(tableName);
            for (String row : rows)
                stream.println(row);
            stream.println();
        }
    }

    /**
     * get the header (tag) names
     * @param analysis the analysis node
     * @return a string representing the tag names
     */
    private String getTagValues(Node analysis) {
        StringBuilder buffer = new StringBuilder();
        for (Node s : analysis.getChildren())
            if (s.tag) buffer.append(formatColumn(s.getValue()));
        return buffer.toString();
    }

    /**
     * simple data points describe themselves, and have one child that stores their value and it's description.  Extract the value and
     * convert the list of nodes to a string
     * @param analysis the analysis
     * @return a String representing the values
     */
    private String dataPointNodesToValues(Node analysis) {
        StringBuilder builder = new StringBuilder();
        for (Node n : analysis.getChildren()) {
            if (!n.tag && !n.table) {
                if (n.getChildren().size() > 1) throw new IllegalStateException("Simple data points shouldn't have more than one value");
                if (n.getChildren().size() == 1)
                    builder.append(formatColumn(n.getChildren().iterator().next().getValue()));
            }
        }
        return builder.toString();
    }

    /**
     * extract the header string from the base analysis node
     */
    private String extractHeaderString(Node analysisNode) {
        StringBuilder buffer = new StringBuilder();
        // first get the tags
        getTagNames(analysisNode, buffer);
        if (!getColumnNames(analysisNode, buffer))
            return null;

        return buffer.toString();
    }

    /**
     * get the column names from the analysis node
     * @param analysisNode the node
     * @param buffer the buffer to append to
     * @return true if there was data fields to output, false if we dont add data to the column header list
     */
    private boolean getColumnNames(Node analysisNode, StringBuilder buffer) {
        // now get the simple data points
        boolean addedValue = false;
        for (Node n : analysisNode.getChildren())
            if (!n.tag && !n.table) {
                addedValue = true;
                buffer.append(formatColumn(n.getValue()));
            }
        return addedValue;
    }

    /**
     * get the tags names from an analysis node
     * @param analysisNode the node
     * @param buffer the StringBuilder to append to
     */
    private void getTagNames(Node analysisNode, StringBuilder buffer) {
        for (Node n : analysisNode.getChildren())
            if (n.tag) buffer.append(formatColumn(n.getName()));
    }

    /**
     * create a correct-length divider string
     * @param length the length for the divider
     * @return a string with the divider text of length "length"
     */
    private String niceDivider(int length) {
        if (!displayDashedLineBreaks()) return "";
        StringBuilder builder = new StringBuilder();
        for (int x = 0; x < length; x++) builder.append("-");
        return builder.toString();
    }

    /**
     * format the string according to our internal rules
     * @param str the string to format
     * @return a string, properly formatted
     */
    public abstract String formatColumn(String str);

    /**
     * does the output format want to display line breaks (dotted lines)?
     * @return true if the format uses them
     */
    public abstract boolean displayDashedLineBreaks();
}

package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;

import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class RFormat
 *
 * a format for outputting R data - experimental
 */
public class RFormat implements ReportFormat {
    private Map<String, List<Node>> analyses = new HashMap<String, List<Node>>();
    private PrintWriter stream;

    // some reused R table strings
    private String tempTableName = "temp";
    private String toTableString = " <- as.table(temp)";
    private String tableString = tempTableName + " <- matrix(c(";
    private String rowNamesString = "rownames(" + tempTableName + ") <- c(";
    private String colNamesString = "colnames(" + tempTableName + ") <- c(";

    /**
     * write the analyses out
     * @param writeTo the writer to write to
     * @param baseNode the base node, containing the analyses
     */
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
     *
     * @param baseNode the root node
     */
    private void getAnalyses(Node baseNode) {
        for (Node n : baseNode.getChildren())
            if (!n.tag && n.getComplex()) {
                if (!analyses.containsKey(n.getValue()))
                    analyses.put(n.getValue(), new ArrayList<Node>());
                analyses.get(n.getValue()).add(n);
            }
    }

    /**
     * write the analysis nodes out, only outputting the simple data points (non-table data)
     *
     * @param nodes a list of nodes, of the same analysis type
     */
    private void writeAnalysis(List<Node> nodes) {
        if (nodes == null || nodes.size() < 1 || !nodes.get(0).getName().equals("analysis")) return;
        Node forTitle = nodes.get(0);

        String header = extractHeaderString(forTitle);
        if (header == null) return; // a null here indicates we don't have any unique columns to display

        StringBuilder rowNameBuilder = new StringBuilder();
        StringBuilder rowValueBuilder = new StringBuilder();
        for (Node analysis : nodes) {
            String dataString = dataPointNodesToValues(analysis);
            if (dataString.length() > 0 && !dataString.equals("<null>,")) {
                rowNameBuilder.append("\""+getTagValues(analysis).substring(0,getTagValues(analysis).length() - 1)+"\",");
                rowValueBuilder.append(dataString);
            }            
        }
        stream.print(tableString + rowValueBuilder.toString().substring(0, rowValueBuilder.toString().length() - 1) + "),");
        stream.println("ncol=" + header.split(",").length + ")");
        stream.println(rowNamesString + escapeString(rowNameBuilder.toString().substring(0, rowNameBuilder.toString().length() - 1)) + ")");
        stream.println(colNamesString + escapeString(header.substring(0, header.length() - 1)) + ")");

        stream.println(forTitle.getValue().replace(" ", "_").replace("/","_to_") + toTableString);
        stream.println();
    }

    /**
     * output the tables: look at list of analysis nodes (all from the same analysis) and output the table
     *
     * @param nodes the list of analysis nodes (of the same underlying type)
     */
    void outputTables(List<Node> nodes) {
        Map<String, String> tableRowValues = new HashMap<String, String>();
        Map<String, String> tableRowNames = new HashMap<String, String>();
        Map<String, String> tableHeaders = new HashMap<String, String>();
        Map<String, Integer> columnSizes = new HashMap<String, Integer>();
        for (Node analysis : nodes)

            for (Node n : analysis.getChildren()) {
                if (n.table) {
                    StringBuilder columnNamesBuilder = new StringBuilder();
                    StringBuilder rowValueBuilder = new StringBuilder();
                    StringBuilder rowNamesBuilder = new StringBuilder();

                    for (Node row : n.getChildren()) {
                        int colSize = 0;
                        rowNamesBuilder.append("\""+getTagValues(analysis));
                        rowNamesBuilder.append(row.getValue()+"\",");
                        for (Node column : row.getChildren()) {
                            columnNamesBuilder.append("\""+column.getValue()+"\",");
                            if (column.getChildren().size() == 1) {
                                String value = column.getChildren().iterator().next().getValue()+",";
                                colSize++;
                                rowValueBuilder.append(value);
                            }
                        }
                        if (!columnSizes.containsKey(n.getValue())) columnSizes.put(n.getValue(),colSize);
                        if (!tableHeaders.containsKey(n.getValue())) {
                            tableHeaders.put(n.getValue(), columnNamesBuilder.toString());
                        }

                    }
                    tableRowValues.put(n.getValue(),rowValueBuilder.toString());
                    tableRowNames.put(n.getValue(),rowNamesBuilder.toString());
                }
            }

        // output the tables
        for (String tableName : tableHeaders.keySet()) {
            stream.print(tableString + tableRowValues.get(tableName).substring(0,tableRowValues.get(tableName).length()-1) + "),");
            stream.println("ncol="+columnSizes.get(tableName)+")");
            stream.println(rowNamesString + escapeString(tableRowNames.get(tableName).substring(0,tableRowNames.get(tableName).length()-1)) + ")");
            stream.println(colNamesString + escapeString(tableHeaders.get(tableName).substring(0,tableHeaders.get(tableName).length()-1)) + ")");

            stream.println(tableName.replace(" ","_") + toTableString);
            stream.println();
        }
    }

    /**
     * get the header (tag) names
     *
     * @param analysis the analysis node
     *
     * @return a string representing the tag names
     */
    private String getTagValues(Node analysis) {
        StringBuilder buffer = new StringBuilder();
        for (Node s : analysis.getChildren())
            if (s.tag) buffer.append(s.getValue()+"_");
        return buffer.toString();
    }

    /**
     * escape any special characters for R
     * @param str the string to check
     * @return the escaped string
     */
    String escapeString(String str) {
        if (str.contains("/")) str = str.replace("/","\\/");
        return str;
    }

    /**
     * convert the list of nodes to a string
     *
     * @param analysis the analysis
     *
     * @return a String representing the values
     */
    private String dataPointNodesToValues(Node analysis) {
        StringBuilder builder = new StringBuilder();
        for (Node n : analysis.getChildren()) {
            if (!n.tag && !n.table) {
                if (n.getChildren().size() > 1)
                    throw new IllegalStateException("Simple data points shouldn't have more than one value");
                if (n.getChildren().size() == 1)
                    builder.append(formatColumn(n.getChildren().iterator().next().getValue()));
            }
        }
        return builder.toString();
    }

    /** extract the header string from the base analysis node */
    private String extractHeaderString(Node analysisNode) {
        StringBuilder buffer = new StringBuilder();
        // first get the tags
        if (!getColumnNames(analysisNode, buffer))
            return null;

        return buffer.toString();
    }

    /**
     * get the column names from the analysis node
     *
     * @param analysisNode the node
     * @param buffer       the buffer to append to
     *
     * @return true if there was data fields to output, false if we don't add data to the column header list
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
     * format a column string
     * @param row the column
     * @return the reformatted string
     */
    public String formatColumn(String row) { return "\"" + row + "\",";}
}

package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.exceptions.UserError;

import java.io.*;
import java.util.*;

/**
 * an abstract class to share the basics of a table based format; many methods
 * overlap in different output types.
 */
public abstract class TableBasedFormat implements ReportFormat {
    private Map<String, List<Node>> analyses = new HashMap<String, List<Node>>();
    private PrintWriter stream;
    private File baseLocation;

    /**
     * write the base node to the specified writer
     * @param writeTo the file base to write to
     * @param baseNode the root node
     */
    @Override
    public void write(File writeTo, Node baseNode) {
        baseLocation = writeTo;

        // if there is only a single output file, create it
        if (!splitFilesByAnalysis()) newStream("");

        traverseAnalysisNodes(baseNode);
    }

    /**
     * write the base node to the specified writer
     * @param writeLocation the writer to write to
     * @param baseNode the root node
     */
    public void write(Writer writeLocation, Node baseNode) {
        if (splitFilesByAnalysis()) throw new UserError.CommandLineError("Unable to write output report, we require a file input for multi-file formats");
        // if there is only a single output file, create it
        stream = new PrintWriter(writeLocation);
        traverseAnalysisNodes(baseNode);
        stream.flush();
    }

    /**
     * traverse the analysis nodes, outputting to our stream
     * @param baseNode the base (root) node, with analysis nodes as children
     */
    private void traverseAnalysisNodes(Node baseNode) {
        getAnalyses(baseNode);
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
        newStream(forTitle.getValue());
        stream.println(headerIndicator() + "Analysis Name:         \t" + forTitle.getValue());
        stream.println(headerIndicator() + "Analysis Description:  \t" + forTitle.getDescription());
        if (addReadabilityMarks()) stream.println();

        String header = extractHeaderString(forTitle);
        if (header == null) return; // a null here indicates we don't have any unique columns to display
        stream.println(trimLastChar(header));
        if (addReadabilityMarks()) stream.println(niceDivider(header.length()));

        for (Node analysis : nodes) {
            String dataString = dataPointNodesToValues(analysis);
            if (dataString.length() > 0 && !dataString.equals("<null>")) {
                stream.print(getTagValues(analysis));
                stream.println(trimLastChar(dataString));
            }
        }
        if (addReadabilityMarks()) stream.println();
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
            newStream(tableName);
            stream.println(headerIndicator() + "Table Name : " + tableName);
            stream.println(trimLastChar(tableHeaders.get(tableName)));
            if (addReadabilityMarks()) stream.println(niceDivider(tableHeaders.get(tableName).length()));
            List<String> rows = tableRows.get(tableName);
            for (String row : rows)
                stream.println(trimLastChar(row));            
            if (addReadabilityMarks()) stream.println();
        }
    }

    public String trimLastChar(String toTrim) {
        return toTrim.substring(0,toTrim.length()-1);
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
     * this function checks whether we need to create a new stream for the specified analysis
     */
    public void newStream(String analysisOrTableName) {
        String name = analysisOrTableName.replaceAll("\\s+","_").replaceAll("\\/","_slash_");
        File file = new File(this.baseLocation + "." + name + this.extension());
        if (stream == null || splitFilesByAnalysis()) {
            if (stream != null) stream.close();
            try {
                stream = new PrintWriter(file);
            } catch (FileNotFoundException e) {
                throw new UserError.CouldNotCreateOutputFile(file, e);
            }
        }
    }

    /**
     * return the valid outputs we support
     * @return
     */
    public EnumSet<AcceptableOutputType> getAcceptableOutputTypes() {
        EnumSet<AcceptableOutputType> set =  EnumSet.of(AcceptableOutputType.FILE); // always acceptable
        if (!splitFilesByAnalysis()) set.add(AcceptableOutputType.STREAM);
        return set;
    }

    /**
     * create a correct-length divider string
     * @param length the length for the divider
     * @return a string with the divider text of length "length"
     */
    private String niceDivider(int length) {
        StringBuilder builder = new StringBuilder();
        for (int x = 0; x < length; x++) builder.append("-");
        return builder.toString();
    }

    /**
     * close the output file, if open
     */
    public void close() {
        if (stream != null) stream.close();
    }
    
    /**
     * format the string according to our internal rules
     * @param str the string to format
     * @return a string, properly formatted
     */
    public abstract String formatColumn(String str);

    /**
     * should we add readability marks?
     * @return true if we should (line breaks, etc)
     */
    public abstract boolean addReadabilityMarks();


    /**
     * a string to prepend for header lines
     * @return a string, blank if no string to be appended
     */
    public abstract String headerIndicator();

    /**
     * should we split the separate files by analysis
     * @return
     */
    public abstract boolean splitFilesByAnalysis();

    /**
     * what extension do we want our files to have
     * @return a string of the extension
     */
    public abstract String extension();
}

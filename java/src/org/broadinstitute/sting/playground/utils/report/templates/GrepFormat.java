package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.EnumSet;


/**
 * 
 * @author aaron 
 * 
 * Class GrepFormat
 *
 * implements the grep output format
 */
public class GrepFormat implements ReportFormat {
    private PrintWriter stream;

    /**
     * write out to the writer, given the root node
     * @param baseFile the file to write to
     * @param baseNode the root node
     */
    @Override
    public void write(File baseFile, Node baseNode) {
        try {
            stream = new PrintWriter(baseFile);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(baseFile, e);
        }
        privateWrite(baseNode);
    }

    /**
     * write out to the writer, given the root node
     *
     * @param baseFile the file to write to
     * @param baseNode the root node
     */
    @Override
    public void write(Writer baseFile, Node baseNode) {
        stream = new PrintWriter(baseFile);
        privateWrite(baseNode);
    }

    /**
     * write out the node
     * @param baseNode the base (root) node
     */
    private void privateWrite(Node baseNode) {
        for (Node analysis : baseNode.getChildren()) {
            StringBuilder builder = new StringBuilder();
            boolean first = true;
            for (Node tag : analysis.getChildren()) {
                if (first) first = false;
                else if (tag.tag) {
                    builder.append(".");
                }
                if ( tag.tag ) builder.append("["+tag.getName() + "=" + tag.getValue()+"]");
            }
            recursiveTraverse(analysis,builder.toString());
        }
    }

    /**
     * recursively get the data.  If we hit a final node, output the node plus the trailing text
     * @param n the node we're looking at
     * @param value the previous text we've seen
     */
    public void recursiveTraverse(Node n, String value) {
        if (n.tag) return;
        if (n.getChildren().size() < 1) {
            stream.println(value + " " + n.getValue());
        }
        else {
            String nString = n.getName() + "=" +n.getValue();
            for (Node child : n.getChildren())
                recursiveTraverse(child,value + ".[" + nString + "]");
        }
    }

    @Override
    public void close() {
        stream.close();
    }

    /**
     * return the valid outputs we support
     * @return
     */
    public EnumSet<AcceptableOutputType> getAcceptableOutputTypes() {
        EnumSet<AcceptableOutputType> set =  EnumSet.of(AcceptableOutputType.FILE); // always acceptable
        set.add(AcceptableOutputType.STREAM);
        return set;
    }
}

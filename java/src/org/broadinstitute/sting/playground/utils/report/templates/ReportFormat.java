package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;

import java.io.File;
import java.io.Writer;
import java.util.EnumSet;

/**
 * @author aaron
 *         <p/>
 *         Interface ReportFormat
 *         <p/>
 *         The basics of a report formatter
 */
public interface ReportFormat {
    public enum AcceptableOutputType { STREAM, FILE };
    public EnumSet<AcceptableOutputType> getAcceptableOutputTypes();
    public void write(File fileLocation, Node baseNode);
    public void write(Writer writeLocation, Node baseNode);
    public void close();
}

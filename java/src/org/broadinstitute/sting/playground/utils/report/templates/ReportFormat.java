package org.broadinstitute.sting.playground.utils.report.templates;

import org.broadinstitute.sting.playground.utils.report.utils.Node;

import java.io.File;
import java.io.Writer;

/**
 * @author aaron
 *         <p/>
 *         Interface ReportFormat
 *         <p/>
 *         The basics of a report formatter
 */
public interface ReportFormat {
    public void write(Writer baseFile, Node baseNode);
}

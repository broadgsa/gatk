/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.diffengine;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportColumn;
import org.broadinstitute.sting.gatk.report.GATKReportTable;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;


/**
 * Class implementing diffnode reader for GATKReports
 */
public class GATKReportDiffableReader implements DiffableReader {
    @Override
    public String getName() { return "GATKReport"; }

    @Override
    public DiffElement readFromFile(File file, int maxElementsToRead) {
        DiffNode root = DiffNode.rooted(file.getName());
        try {
            // one line reads the whole thing into memory
            GATKReport report = new GATKReport(file);

            for (GATKReportTable table : report.getTables() ) {
                root.add(tableToNode(table, root));
            }

            return root.getBinding();
        } catch ( Exception e ) {
            return null;
        }
    }

    private DiffNode tableToNode(GATKReportTable table, DiffNode root) {
        DiffNode tableRoot = DiffNode.empty(table.getTableName(), root);

        tableRoot.add("Description", table.getTableDescription());
        tableRoot.add("NumberOfRows", table.getNumRows());
        tableRoot.add("Version", table.getVersion());

        for ( GATKReportColumn column : table.getColumns().values() ) {
            DiffNode columnRoot = DiffNode.empty(column.getColumnName(), tableRoot);

            columnRoot.add("Width", column.getColumnWidth());
            columnRoot.add("Displayable", column.isDisplayable());

            int n = 1;
            for ( Object elt : column.values() ) {
                String name = column.getColumnName() + n++;
                columnRoot.add(name, elt.toString());
            }

            tableRoot.add(columnRoot);
        }

        return tableRoot;
    }

    @Override
    public boolean canRead(File file) {
        try {
            final String HEADER = GATKReport.GATKREPORT_HEADER_PREFIX;
            char[] buff = new char[HEADER.length()];
            new FileReader(file).read(buff, 0, HEADER.length());
            String firstLine = new String(buff);
            return firstLine.startsWith(HEADER);
        } catch ( IOException e ) {
            return false;
        }
    }
}

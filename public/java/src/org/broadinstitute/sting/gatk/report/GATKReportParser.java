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

package org.broadinstitute.sting.gatk.report;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class GATKReportParser {
    private List<GATKReportTableParser> tables = new ArrayList<GATKReportTableParser>();

    public void parse(File file) throws IOException {
        InputStream stream = FileUtils.openInputStream(file);
        try {
            parse(stream);
        } finally {
            IOUtils.closeQuietly(stream);
        }
    }

    public void parse(InputStream input) throws IOException {
        GATKReportTableParser table = null;

        for (String line: new XReadLines(input)) {
            if (line.startsWith("##:GATKReport.v0.1 ")) {
                table = newTableParser(line);
                tables.add(table);
                table.parse(line);
            } else if (table != null) {
                if (line.trim().length() == 0)
                    table = null;
                else
                    table.parse(line);
            }
        }
    }

    public String getValue(String tableName, String[] key, String column) {
        for (GATKReportTableParser table: tables)
            if (table.getTableName().equals(tableName))
                return table.getValue(key, column);
        return null;
    }

    public String getValue(String tableName, String key, String column) {
        for (GATKReportTableParser table: tables)
            if (table.getTableName().equals(tableName))
                return table.getValue(key, column);
        return null;
    }

    private GATKReportTableParser newTableParser(String header) {
        return new GATKReportTableParser();
    }
}

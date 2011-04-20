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

import org.apache.commons.lang.StringUtils;

import java.util.*;

public class GATKReportTableParser {
    private int lineNum = 0;
    private String[] descriptions;
    private Map<String, Integer> headers = new HashMap<String, Integer>();
    private List<String[]> values = new ArrayList<String[]>();

    public void parse(String line) {
        lineNum++;
        switch (lineNum) {
            case 1:
                descriptions = parseLine(line);
            case 2:
                String[] columnHeaders = parseLine(line);
                for (int i = 0; i < columnHeaders.length; i++)
                    headers.put(columnHeaders[i], i);
            default:
                values.add(parseLine(line));
        }
    }

    public String getTableName() {
        return descriptions[1];
    }

    public String getValue(String[] key, String column) {
        if (!headers.containsKey(column))
            return null;
        for (String[] row: values)
            if (Arrays.equals(key, Arrays.copyOfRange(row, 1, key.length + 1)))
                return row[headers.get(column)];
        return null;
    }

    public String getValue(String key, String column) {
        return getValue(key.split("\\."), column);
    }

    private String generateKey(String[] row, int i) {
        return StringUtils.join(row, ".", 0, i);
    }

    private String[] parseLine(String line) {
        return line.split(" +");
    }
}

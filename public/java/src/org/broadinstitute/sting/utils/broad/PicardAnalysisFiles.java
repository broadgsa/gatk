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

package org.broadinstitute.sting.utils.broad;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class PicardAnalysisFiles {
    private static final String REFERENCE_SEQUENCE_HEADER = "REFERENCE_SEQUENCE";
    private static final String TARGET_INTERVALS_HEADER = "TARGET_INTERVALS";
    private static final String BAIT_INTERVALS_HEADER = "BAIT_INTERVALS";
    private static final String[] ANALYSIS_HEADERS = {REFERENCE_SEQUENCE_HEADER, TARGET_INTERVALS_HEADER, BAIT_INTERVALS_HEADER};
    private static final String ANALYSIS_FILES = "analysis_files.txt";

    private String path;
    private Map<String,Set<String>> headerValues = new HashMap<String,Set<String>>();

    public PicardAnalysisFiles(String project, String sample) throws FileNotFoundException {
        this(PicardAggregationUtils.getSampleDir(project, sample) + ANALYSIS_FILES);
    }

    public PicardAnalysisFiles(String project, String sample, int version) throws FileNotFoundException {
        this(PicardAggregationUtils.getSampleDir(project, sample, version) + ANALYSIS_FILES);
    }

    public PicardAnalysisFiles(String path) throws FileNotFoundException {
        this.path = path;
        HashMap<String,Integer> headerIndexes = null;
        for (String line: new XReadLines(new File(path))) {
            if (line.startsWith("#"))
                continue;
            String[] values = line.split("\t");
            if (headerIndexes == null) {
                headerIndexes = new HashMap<String,Integer>();
                for (String header: ANALYSIS_HEADERS) {
                    headerIndexes.put(header, ArrayUtils.indexOf(values, header));
                    headerValues.put(header, new HashSet<String>());
                }
            } else {
                for (String header: ANALYSIS_HEADERS) {
                    int index = headerIndexes.get(header);
                    if (values.length <= index)
                        throw new StingException(String.format("Unable to parse line in %s: %n%s", path, line));
                    String value = values[index];
                    headerValues.get(header).add(value);
                }
            }
        }
    }

    public String getPath() {
        return path;
    }

    public String getReferenceSequence() {
        return getSingle(REFERENCE_SEQUENCE_HEADER);
    }

    public String getTargetIntervals() {
        return getSingle(TARGET_INTERVALS_HEADER);
    }

    public String getBaitIntervals() {
        return getSingle(BAIT_INTERVALS_HEADER);
    }

    private String getSingle(String header) {
        Set<String> values = headerValues.get(header);
        if (values.size() > 1) {
            throw new UnsupportedOperationException(path + " contains more than one value for " + header + ": " + values);
        } else if (values.size() == 0) {
            return null;
        } else {
            String value = values.iterator().next();
            return "null".equals(value) ? null : value;
        }
    }
}

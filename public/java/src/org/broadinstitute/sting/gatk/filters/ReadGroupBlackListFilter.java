/*
 * Copyright (c) 2010 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.Map.Entry;

/**
 * Removes records matching the read group tag and exact match string.
 * For example, this filter value:
 *   PU:1000G-mpimg-080821-1_1
 * would filter out a read with the read group PU:1000G-mpimg-080821-1_1
 */
public class ReadGroupBlackListFilter extends ReadFilter {
    private Set<Entry<String, Collection<String>>> filterEntries;

    public ReadGroupBlackListFilter(List<String> blackLists) {
        Map<String, Collection<String>> filters = new TreeMap<String, Collection<String>>();
        for (String blackList : blackLists)
            addFilter(filters, blackList, null, 0);
        this.filterEntries = filters.entrySet();
    }

    public boolean filterOut(SAMRecord samRecord) {
        for (Entry<String, Collection<String>> filterEntry : filterEntries) {
            String attributeType = filterEntry.getKey();

            SAMReadGroupRecord samReadGroupRecord = samRecord.getReadGroup();
            if (samReadGroupRecord != null) {
                Object attribute;
                if ("ID".equals(attributeType) || "RG".equals(attributeType))
                    attribute = samReadGroupRecord.getId();
                else
                    attribute = samReadGroupRecord.getAttribute(attributeType);
                if (attribute != null && filterEntry.getValue().contains(attribute))
                    return true;
            }
        }

        return false;
    }

    private void addFilter(Map<String, Collection<String>> filters, String filter, File parentFile, int parentLineNum) {
        if (filter.toLowerCase().endsWith(".list") || filter.toLowerCase().endsWith(".txt")) {
            File file = new File(filter);
            try {
                int lineNum = 0;
                XReadLines lines = new XReadLines(file);
                for (String line : lines) {
                    lineNum++;

                    if (line.trim().length() == 0)
                        continue;

                    if (line.startsWith("#"))
                        continue;

                    addFilter(filters, line, file, lineNum);
                }
            } catch (FileNotFoundException e) {
                String message = "Error loading black list: " + file.getAbsolutePath();
                if (parentFile != null) {
                    message += ", " + parentFile.getAbsolutePath() + ":" + parentLineNum;
                }
                throw new UserException(message);
            }
        } else {
            String[] filterEntry = filter.split(":", 2);

            String message = null;
            if (filterEntry.length != 2) {
                message = "Invalid read group filter: " + filter;
            } else if (filterEntry[0].length() != 2) {
                message = "Tag is not two characters: " + filter;
            }

            if (message != null) {
                if (parentFile != null) {
                    message += ", " + parentFile.getAbsolutePath() + ":" + parentLineNum;
                }
                message += ", format is <TAG>:<SUBSTRING>";
                throw new UserException(message);
            }

            if (!filters.containsKey(filterEntry[0]))
                filters.put(filterEntry[0], new TreeSet<String>());
            filters.get(filterEntry[0]).add(filterEntry[1]);
        }
    }
}

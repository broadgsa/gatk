/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import org.apache.commons.collections.CollectionUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.Gatherer;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * User: carneiro
 * Date: 3/29/11
 */


public class BQSRGatherer extends Gatherer  {

    private static final Logger logger = Logger.getLogger(BQSRGatherer.class);
    private static final String EMPTY_INPUT_LIST = "list of inputs files is empty or there is no usable data in any input file";
    private static final String MISSING_OUTPUT_FILE = "missing output file name";
    private static final String MISSING_READ_GROUPS = "Missing read group(s)";

    @Override
    public void gather(final List<File> inputs, final File output) {
        final PrintStream outputFile;
        try {
            outputFile = new PrintStream(output);
        } catch(FileNotFoundException e) {
            throw new UserException.MissingArgument("output", MISSING_OUTPUT_FILE);
        }
        final GATKReport report = gatherReport(inputs);
        report.print(outputFile);
    }

    /**
     * Gathers the input recalibration reports into a single report.
     *
     * @param inputs Input recalibration GATK reports
     * @return gathered recalibration GATK report
     */
    public static GATKReport gatherReport(final List<File> inputs) {
        final SortedSet<String> allReadGroups = new TreeSet<String>();
        final LinkedHashMap<File, Set<String>> inputReadGroups = new LinkedHashMap<File, Set<String>>();

        // Get the read groups from each input report
        for (final File input : inputs) {
            final Set<String> readGroups = RecalibrationReport.getReadGroups(input);
            inputReadGroups.put(input, readGroups);
            allReadGroups.addAll(readGroups);
        }

        // Log the read groups that are missing from specific inputs
        for (Map.Entry<File, Set<String>> entry: inputReadGroups.entrySet()) {
            final File input = entry.getKey();
            final Set<String> readGroups = entry.getValue();
            if (allReadGroups.size() != readGroups.size()) {
                // Since this is not completely unexpected, more than debug, but less than a proper warning.
                logger.info(MISSING_READ_GROUPS + ": " + input.getAbsolutePath());
                for (final Object readGroup: CollectionUtils.subtract(allReadGroups, readGroups)) {
                    logger.info("  " + readGroup);
                }
            }
        }

        RecalibrationReport generalReport = null;
        for (File input : inputs) {
            final RecalibrationReport inputReport = new RecalibrationReport(input, allReadGroups);
            if( inputReport.isEmpty() ) { continue; }

            if (generalReport == null)
                generalReport = inputReport;
            else
                generalReport.combine(inputReport);
        }
        if (generalReport == null)
            throw new ReviewedGATKException(EMPTY_INPUT_LIST);

        generalReport.calculateQuantizedQualities();

        return generalReport.createGATKReport();
    }
}

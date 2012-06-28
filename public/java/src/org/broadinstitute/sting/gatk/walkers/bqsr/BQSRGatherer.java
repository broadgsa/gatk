/*
 * Copyright (c) 2011 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.commandline.Gatherer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * User: carneiro
 * Date: 3/29/11
 */


public class BQSRGatherer extends Gatherer  {
    
    private static final String EMPTY_INPUT_LIST = "list of inputs files is empty";
    private static final String MISSING_OUTPUT_FILE = "missing output file name";

    @Override
    public void gather(List<File> inputs, File output) {
        RecalibrationReport generalReport = null;
        PrintStream outputFile;
        try {
            outputFile = new PrintStream(output);
        } catch(FileNotFoundException e) {
            throw new UserException.MissingArgument("output", MISSING_OUTPUT_FILE);
        }

        for (File input : inputs) {
            RecalibrationReport inputReport = new RecalibrationReport(input);
            if (generalReport == null)
                generalReport = inputReport;
            else
                generalReport.combine(inputReport);
        }
        if (generalReport == null)
            throw new ReviewedStingException(EMPTY_INPUT_LIST);

        generalReport.calculateEmpiricalAndQuantizedQualities();

        RecalibrationArgumentCollection RAC = generalReport.getRAC();
        if (RAC.recalibrationReport != null && !RAC.NO_PLOTS) {
            File recal_out = new File(output.getName() + ".original");
            RecalibrationReport originalReport = new RecalibrationReport(RAC.recalibrationReport);
            // TODO -- fix me
            //RecalDataManager.generateRecalibrationPlot(recal_out, originalReport.getKeysAndTablesMap(), generalReport.getKeysAndTablesMap(), RAC.KEEP_INTERMEDIATE_FILES);
        }
        else if (!RAC.NO_PLOTS) {
            File recal_out = new File(output.getName() + ".recal");
            // TODO -- fix me
            //RecalDataManager.generateRecalibrationPlot(recal_out, generalReport.getKeysAndTablesMap(), RAC.KEEP_INTERMEDIATE_FILES);
        }

        generalReport.output(outputFile);
    }
}

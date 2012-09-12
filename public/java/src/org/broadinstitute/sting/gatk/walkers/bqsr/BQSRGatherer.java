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
import org.broadinstitute.sting.utils.recalibration.RecalUtils;
import org.broadinstitute.sting.utils.recalibration.RecalibrationReport;

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
        final PrintStream outputFile;
        try {
            outputFile = new PrintStream(output);
        } catch(FileNotFoundException e) {
            throw new UserException.MissingArgument("output", MISSING_OUTPUT_FILE);
        }

        RecalibrationReport generalReport = null;
        for (File input : inputs) {
            final RecalibrationReport inputReport = new RecalibrationReport(input);
            if (generalReport == null)
                generalReport = inputReport;
            else
                generalReport.combine(inputReport);
        }
        if (generalReport == null)
            throw new ReviewedStingException(EMPTY_INPUT_LIST);

        generalReport.calculateQuantizedQualities();

        RecalibrationArgumentCollection RAC = generalReport.getRAC();
        if ( RAC.RECAL_PDF_FILE != null ) {
            RAC.RECAL_TABLE_FILE = output;
            if ( RAC.existingRecalibrationReport != null ) {
                final RecalibrationReport originalReport = new RecalibrationReport(RAC.existingRecalibrationReport);
                RecalUtils.generateRecalibrationPlot(RAC, originalReport.getRecalibrationTables(), generalReport.getRecalibrationTables(), generalReport.getCovariates());
            }
            else {
                RecalUtils.generateRecalibrationPlot(RAC, generalReport.getRecalibrationTables(), generalReport.getCovariates());
            }
        }

        generalReport.output(outputFile);
    }
}

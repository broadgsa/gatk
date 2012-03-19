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
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatumOptimized;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * User: carneiro
 * Date: 3/29/11
 */


public class BQSRGatherer extends Gatherer  {

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private static final String EOF_MARKER = "EOF";

    private HashMap<String, RecalDatumOptimized> dataMap = new HashMap<String, RecalDatumOptimized>();


    private void addCSVData (String line) {
        String[] covariates = line.split(",");
        String key = "";
        RecalDatumOptimized values;

        for (int i = 0; i < covariates.length-3; i++)
            key += covariates[i] + ",";

        if (covariates.length < 3)
            throw new ReviewedStingException("Line only has 1 covariate : " + line);

        values = new RecalDatumOptimized(Long.parseLong(covariates[covariates.length - 3]), Long.parseLong(covariates[covariates.length - 2]));

        RecalDatumOptimized currentValues = dataMap.get(key);
        if (currentValues == null)
            dataMap.put(key, values);
        else
            currentValues.increment(values);

    }

    @Override
    public void gather(List<File> inputs, File output) {
        PrintStream o;
        try {
            o = new PrintStream(output);
        } catch ( FileNotFoundException e) {
            throw new UserException("File to be output by CountCovariates Gather function was not found");
        }

        boolean sawEOF = false;
        boolean printedHeader = false;

        // Read input files
        for ( File RECAL_FILE : inputs) {
            try {
                for ( String line : new XReadLines(RECAL_FILE) ) {
                    if ( EOF_MARKER.equals(line) ) {
                        sawEOF = true;                      // sanity check
                        break;
                    }

                    else if(line.startsWith("#"))  {
                        if (!printedHeader)
                            o.println(line);
                    }

                    else                                    // Found a line of data
                        addCSVData(line);                   // Parse the line and add the data to the HashMap
                }

            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotReadInputFile(RECAL_FILE, "Can not find input file", e);
            }

            if ( !sawEOF ) {
                final String errorMessage = "No EOF marker was present in the recal covariates table; this could mean that the file is corrupted!";
                throw new UserException.MalformedFile(RECAL_FILE, errorMessage);
            }
            printedHeader = true;
        }

        // Write output file from dataMap
        for(Map.Entry<String, RecalDatumOptimized> entry : dataMap.entrySet())
            o.println(entry.getKey() + entry.getValue().outputToCSV());
        o.println("EOF");
    }
}

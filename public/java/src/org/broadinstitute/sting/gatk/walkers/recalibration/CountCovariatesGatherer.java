package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.commandline.Gatherer;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 3/29/11
 * Time: 3:54 PM
 * To change this template use File | Settings | File Templates.
 */


public class CountCovariatesGatherer extends Gatherer  {

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private static final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private static final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");
    private static final String EOF_MARKER = "EOF";

    private HashMap<String, RecalDatumOptimized> dataMap;


    private void addCSVData (String line) {
        String[] covariates = line.split(",");
        String key = "";
        RecalDatumOptimized values;

        for (int i = 0; i < covariates.length-3; i++) {
            key += covariates[i] + ",";
        }

        values = new RecalDatumOptimized(Integer.parseInt(covariates[covariates.length-3]),
                                         Integer.parseInt(covariates[covariates.length-2]));

        if (dataMap.get(key) != null) {
            RecalDatumOptimized currentValues = dataMap.get(key);
            values.increment(currentValues);
        }

        dataMap.put(key, values);
    }

    @Override
    public void gather(List<File> inputs, File output) {
        dataMap = new HashMap<String, RecalDatumOptimized>();
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
                        sawEOF = true;    // sanity check
                    }
                    else if(COMMENT_PATTERN.matcher(line).matches()) {
                        ;                 // It doesn't make any sense to print intermediate comments, unless we merge them somehow (would require strict definition for the header)
                    }
                    else if (COVARIATE_PATTERN.matcher(line).matches()) {
                        if (!printedHeader)
                            o.println(line);
                    }
                    else {                // Found a line of data
                        addCSVData(line); // Parse the line and add the data to the HashMap
                    }
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
        for(String key : dataMap.keySet()) {
            RecalDatumOptimized values = dataMap.get(key);
            String v = values.getNumObservations() + "," + values.getNumMismatches() + "," + values.empiricalQualByte();
            o.println(key + v);
        }
        o.println("EOF");
    }
}

/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.analyzecovariates;

import org.broadinstitute.sting.gatk.walkers.recalibration.*;
import org.broadinstitute.sting.utils.classloader.PackageUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Dec 1, 2009
 *
 * Create collapsed versions of the recal csv file and call R scripts to plot residual error versus the various covariates.
 */

class AnalyzeCovariatesCLP extends CommandLineProgram {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////

    @Argument(fullName = "recal_file", shortName = "recalFile", doc = "The input recal csv file to analyze", required = false)
    private String RECAL_FILE = "output.recal_data.csv";
    @Argument(fullName = "output_dir", shortName = "outputDir", doc = "The directory in which to output all the plots and intermediate data files", required = false)
    private String OUTPUT_DIR = "analyzeCovariates/";
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is probably /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";
    @Argument(fullName = "ignoreQ", shortName = "ignoreQ", doc = "Ignore bases with reported quality less than this number.", required = false)
    private int IGNORE_QSCORES_LESS_THAN = 5;
    @Argument(fullName = "numRG", shortName = "numRG", doc = "Only process N read groups. Default value: -1 (process all read groups)", required = false)            
    private int NUM_READ_GROUPS_TO_PROCESS = -1; // -1 means process all read groups
    @Argument(fullName="max_quality_score", shortName="maxQ", required = false, doc="The integer value at which to cap the quality scores, default is 40")
    private int MAX_QUALITY_SCORE = 40;
    @Argument(fullName="max_histogram_value", shortName="maxHist", required = false, doc="If supplied, this value will be the max value of the histogram plots")
    private int MAX_HISTOGRAM_VALUE = 0;


    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private AnalysisDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private ArrayList<Covariate> requestedCovariates; // List of covariates to be used in this calculation
    private final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private final Pattern OLD_RECALIBRATOR_HEADER = Pattern.compile("^rg,.*");
    private final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");

    protected int execute() {

        // create the output directory where all the data tables and plots will go
        try {
            Process p = Runtime.getRuntime().exec("mkdir " + OUTPUT_DIR);
        } catch (IOException e) {
            System.out.println("Couldn't create directory: " + OUTPUT_DIR);
            System.out.println("User is responsible for making sure the output directory exists.");
        }
        if( !OUTPUT_DIR.endsWith("/") ) { OUTPUT_DIR = OUTPUT_DIR + "/"; }
        if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }

        // initialize all the data from the csv file and allocate the list of covariates
        System.out.println("Reading in input csv file...");
        initializeData();
        System.out.println("...Done!");

        // output data tables for Rscript to read in
        System.out.println("Writing out intermediate tables for R...");
        writeDataTables();
        System.out.println("...Done!");

        // perform the analysis using Rscript and output the plots
        System.out.println("Calling analysis R scripts and writing out figures...");
        callRScripts();
        System.out.println("...Done!");

        return 1;
    }

    private void initializeData() {

        // Get a list of all available covariates
        List<Class<? extends Covariate>> classes = PackageUtils.getClassesImplementingInterface(Covariate.class);

        int lineNumber = 0;
        boolean foundAllCovariates = false;

        // Read in the covariates that were used from the input file
        requestedCovariates = new ArrayList<Covariate>();

        try {
            for ( String line : new XReadLines(new File( RECAL_FILE )) ) {
                lineNumber++;
                if( COMMENT_PATTERN.matcher(line).matches() || OLD_RECALIBRATOR_HEADER.matcher(line).matches())  {
                    ; // Skip over the comment lines, (which start with '#')
                }
                else if( COVARIATE_PATTERN.matcher(line).matches() ) { // The line string is either specifying a covariate or is giving csv data
                    if( foundAllCovariates ) {
                        throw new RuntimeException( "Malformed input recalibration file. Found covariate names intermingled with data in file: " + RECAL_FILE );
                    } else { // Found the covariate list in input file, loop through all of them and instantiate them
                        String[] vals = line.split(",");
                        for( int iii = 0; iii < vals.length - 3; iii++ ) { // There are n-3 covariates. The last three items are nObservations, nMismatch, and Qempirical
                            boolean foundClass = false;
                            for( Class<?> covClass : classes ) {
                                if( (vals[iii] + "Covariate").equalsIgnoreCase( covClass.getSimpleName() ) ) {
                                    foundClass = true;
                                    try {
                                        Covariate covariate = (Covariate)covClass.newInstance();
                                        requestedCovariates.add( covariate );

                                    } catch ( InstantiationException e ) {
                                        throw new RuntimeException( String.format("Can not instantiate covariate class '%s': must be concrete class.", covClass.getSimpleName()) );
                                    } catch ( IllegalAccessException e ) {
                                        throw new RuntimeException( String.format("Can not instantiate covariate class '%s': must have no-arg constructor.", covClass.getSimpleName()) );
                                    }
                                }
                            }

                            if( !foundClass ) {
                                throw new RuntimeException( "Malformed input recalibration file. The requested covariate type (" + (vals[iii] + "Covariate") + ") isn't a valid covariate option." );
                            }
                        }

                    }

                } else { // Found a line of data
                    if( !foundAllCovariates ) {

                        foundAllCovariates = true;

                        // At this point all the covariates should have been found and initialized
                        if( requestedCovariates.size() < 2 ) {
                            throw new RuntimeException( "Malformed input recalibration file. Covariate names can't be found in file: " + RECAL_FILE );
                        }

                        // Initialize any covariate member variables using the shared argument collection
                        for( Covariate cov : requestedCovariates ) {
                            cov.initialize( new RecalibrationArgumentCollection() );
                        }

                        // Initialize the data hashMaps
                        dataManager = new AnalysisDataManager( requestedCovariates.size() );

                    }
                    addCSVData(line); // Parse the line and add the data to the HashMap
                }
            }

        } catch ( FileNotFoundException e ) {
            throw new RuntimeException("Can not find input file: " + RECAL_FILE);
        } catch ( NumberFormatException e ) {
            throw new RuntimeException("Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
        }
    }

    private void addCSVData(String line) {
        String[] vals = line.split(",");

        // Check if the data line is malformed, for example if the read group string contains a comma then it won't be parsed correctly
        if( vals.length != requestedCovariates.size() + 3 ) { // +3 because of nObservations, nMismatch, and Qempirical
            throw new RuntimeException("Malformed input recalibration file. Found data line with too many fields: " + line +
                    " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
        }

        Object[] key = new Object[requestedCovariates.size()];
        Covariate cov;
        int iii;
        for( iii = 0; iii < requestedCovariates.size(); iii++ ) {
            cov = requestedCovariates.get( iii );
            key[iii] = cov.getValue( vals[iii] );
        }
        // Create a new datum using the number of observations, number of mismatches, and reported quality score
        RecalDatum datum = new RecalDatum( Long.parseLong( vals[iii] ), Long.parseLong( vals[iii + 1] ), Double.parseDouble( vals[1] ), 0.0 );
        // Add that datum to all the collapsed tables which will be used in the sequential calculation
        dataManager.addToAllTables( key, datum, IGNORE_QSCORES_LESS_THAN );

    }

    private void writeDataTables() {

        int numReadGroups = 0;

        // for each read group
        for( Object readGroupKey : dataManager.getCollapsedTable(0).data.keySet() ) {

            if(NUM_READ_GROUPS_TO_PROCESS == -1 || ++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS) {
                String readGroup = readGroupKey.toString();
                RecalDatum readGroupDatum = (RecalDatum) dataManager.getCollapsedTable(0).data.get(readGroupKey);
                System.out.print("Writing out data tables for read group: " + readGroup + "\twith " + readGroupDatum.getNumObservations() + " observations"  );
                System.out.println("\tand aggregate residual error = " + String.format("%.3f", readGroupDatum.empiricalQualDouble(0, MAX_QUALITY_SCORE) - readGroupDatum.getEstimatedQReported()));

                // for each covariate
                for( int iii = 1; iii < requestedCovariates.size(); iii++ ) {
                    Covariate cov = requestedCovariates.get(iii);

                    // Create a PrintStream
                    PrintStream output = null;
                    try {
                        output = new PrintStream(new FileOutputStream(OUTPUT_DIR + readGroup + "." + cov.getClass().getSimpleName()+ ".dat"));

                    } catch (FileNotFoundException e) {
                        System.err.println("Can't create file: " + OUTPUT_DIR + readGroup + "." + cov.getClass().getSimpleName()+ ".dat");
                        System.exit(-1);
                    }

                    // Output the header
                    output.println("Covariate\tQreported\tQempirical\tnMismatches\tnBases");

                    for( Object covariateKey : ((Map)dataManager.getCollapsedTable(iii).data.get(readGroupKey)).keySet()) {
                        output.print( covariateKey.toString() + "\t" );                                                     // Covariate
                        RecalDatum thisDatum = (RecalDatum)((Map)dataManager.getCollapsedTable(iii).data.get(readGroupKey)).get(covariateKey);
                        output.print( String.format("%.3f", thisDatum.getEstimatedQReported()) + "\t" );                    // Qreported
                        output.print( String.format("%.3f", thisDatum.empiricalQualDouble(0, MAX_QUALITY_SCORE)) + "\t" );  // Qempirical
                        output.print( thisDatum.getNumMismatches() + "\t" );                                                // nMismatches
                        output.println( thisDatum.getNumObservations() );                                                   // nBases
                    }

                    // Close the PrintStream
                    output.close();
                }
            } else {
                break;
            }

        }
    }

    private void callRScripts() {

        int numReadGroups = 0;
        
        // for each read group
        for( Object readGroupKey : dataManager.getCollapsedTable(0).data.keySet() ) {

            Process p = null;
            if(++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS || NUM_READ_GROUPS_TO_PROCESS == -1) {

                String readGroup = readGroupKey.toString();
                System.out.println("Analyzing read group: " + readGroup);

                // for each covariate
                for( int iii = 1; iii < requestedCovariates.size(); iii++ ) {
                    Covariate cov = requestedCovariates.get(iii);
                    try {
                        if( iii == 1 ) {
                            // Analyze reported quality
                            p = Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_residualError_QualityScoreCovariate.R" + " " +
                                        OUTPUT_DIR + readGroup + "." + cov.getClass().getSimpleName()+ ".dat" + " " +
                                        IGNORE_QSCORES_LESS_THAN + " " + MAX_QUALITY_SCORE + " " + MAX_HISTOGRAM_VALUE); // The third argument is the Q scores that should be turned pink in the plot because they were ignored
                            if(numReadGroups % 3 == 0) { // Don't want to spawn all the RScript jobs too quickly. So wait for this one to finish
                                p.waitFor();
                            }
                        } else { // Analyze all other covariates
                            Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_residualError_OtherCovariate.R" + " " +
                                        OUTPUT_DIR + readGroup + "." + cov.getClass().getSimpleName()+ ".dat" + " " +
                                        cov.getClass().getSimpleName().split("Covariate")[0]); // The third argument is the name of the covariate in order to make the plots look nice
                        }
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                        System.exit(-1);
                    } catch (IOException e) {
                        System.out.println("Fatal Exception: Perhaps RScript jobs are being spawned too quickly? One work around is to process fewer read groups using the -numRG option.");
                        e.printStackTrace();
                        System.exit(-1);
                    }
                }
            } else { // at the maximum number of read groups so break out
                break;
            }
        }
    }
}

public class AnalyzeCovariates {
    public static void main(String args[]) {
        AnalyzeCovariatesCLP clp = new AnalyzeCovariatesCLP();
        CommandLineProgram.start( clp, args );
        System.exit(0);
    }
}

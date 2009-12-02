package org.broadinstitute.sting.analyzecovariates;

import org.broadinstitute.sting.gatk.walkers.recalibration.*;
import org.broadinstitute.sting.utils.PackageUtils;
import org.broadinstitute.sting.utils.xReadLines;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Dec 1, 2009
 */

public class AnalyzeCovariates {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    private static String RECAL_FILE = "output.recal_data.csv";
    private static String OUTPUT_DIR = "analyzeCovariates/";
    private static String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    private static String PATH_TO_RESOURCES = "R/";
    private static int IGNORE_QSCORES_LESS_THAN = 5;
    private static int NUM_READ_GROUPS_TO_PROCESS = -1; // -1 means process all read groups

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private static AnalysisDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private static ArrayList<Covariate> requestedCovariates; // List of covariates to be used in this calculation
    private static final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    private static final Pattern OLD_RECALIBRATOR_HEADER = Pattern.compile("^rg,.*");
    private static final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");

    public static void main(String[] args) {

        // parse command line arguments
        parseArguments( args );
        // create the output directory where all the data tables and plots will go
        try {
            Process p = Runtime.getRuntime().exec("mkdir " + OUTPUT_DIR);
        } catch (IOException e) {
            throw new RuntimeException("Couldn't create directory: " + OUTPUT_DIR);
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

    }

    private static void parseArguments( String[] args ) {
        int iii = 0;
        String arg;

        try {
            while( iii < args.length && args[iii].startsWith("-") ) {
                arg = args[iii++];

                if( arg.equals( "-recalFile" ) ) {
                    RECAL_FILE = args[iii++];
                } else if( arg.equals( "-Rscript" ) ) {
                    PATH_TO_RSCRIPT = args[iii++];
                } else if( arg.equals( "-resources" ) ) {
                    PATH_TO_RESOURCES = args[iii++];
                } else if( arg.equals( "-ignoreQ" ) ) {
                    IGNORE_QSCORES_LESS_THAN = Integer.parseInt( args[iii++] );
                } else if (arg.equals( "-numRG" ) ) {
                    NUM_READ_GROUPS_TO_PROCESS = Integer.parseInt( args[iii++] ); 
                } else if( arg.equals( "-outputDir" ) ) {
                    OUTPUT_DIR = args[iii++];
                } else {
                    iii = -1;
                    break;
                }
            }

            if( iii != args.length ) {
                throw new RuntimeException( "Exception" );
            }
        } catch(Exception e) {
            System.out.println( "Usage: [-option param] \n" );
            System.out.println(" Available options:");
            System.out.println("\t-recalFile <path>\tPath to input recal csv file. Default value: output.recal_data.csv");
            System.out.println("\t-Rscript <path>\t\tPath to your implementation of Rscript. Default value: /broad/tools/apps/R-2.6.0/bin/Rscript");
            System.out.println("\t-resources <path>\tPath to resources folder holding the Sting R scripts. Default value: R/");
            System.out.println("\t-outputDir <path>\tWhere to put the output plots. Default value: analyzeCovariates/");
            System.out.println("\t-ignoreQ <int>\t\tIgnore bases with reported quality less than this number. Default value: 5");
            System.out.println("\t-numRG <int>\t\tOnly process N read groups. Default value: -1 (process all read groups)");
            System.exit(-1);
        }
    }

    private static void initializeData() {

        // Get a list of all available covariates
        List<Class<? extends Covariate>> classes = PackageUtils.getClassesImplementingInterface(Covariate.class);

        int lineNumber = 0;
        boolean foundAllCovariates = false;
        int estimatedCapacity = 1; // Capacity is multiplicitive so this starts at one

        // Read in the covariates that were used from the input file
        requestedCovariates = new ArrayList<Covariate>();

        try {
            for ( String line : new xReadLines(new File( RECAL_FILE )) ) {
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
                                        estimatedCapacity *= covariate.estimatedNumberOfBins();

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

                        // Don't want to crash with out of heap space exception
                        if( estimatedCapacity > 300 * 40 * 200 || estimatedCapacity < 0 ) { // Could be negative if overflowed
                            estimatedCapacity = 300 * 40 * 200;
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

    private static void addCSVData(String line) {
        String[] vals = line.split(",");

        // Check if the data line is malformed, for example if the read group string contains a comma then it won't be parsed correctly
        if( vals.length != requestedCovariates.size() + 3 ) { // +3 because of nObservations, nMismatch, and Qempirical
            throw new RuntimeException("Malformed input recalibration file. Found data line with too many fields: " + line +
                    " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
        }

        ArrayList<Comparable> key = new ArrayList<Comparable>();
        Covariate cov;
        int iii;
        for( iii = 0; iii < requestedCovariates.size(); iii++ ) {
            cov = requestedCovariates.get( iii );
            key.add( cov.getValue( vals[iii] ) );
        }
        // Create a new datum using the number of observations, number of mismatches, and reported quality score
        RecalDatum datum = new RecalDatum( Long.parseLong( vals[iii] ), Long.parseLong( vals[iii + 1] ), Double.parseDouble( vals[1] ) );
        // Add that datum to all the collapsed tables which will be used in the sequential calculation
        dataManager.addToAllTables( key, datum, IGNORE_QSCORES_LESS_THAN );

    }

    private static void writeDataTables() {

        int numReadGroups = 0;

        // for each read group
        NHashMap<RecalDatum> readGroupTable = dataManager.getCollapsedTable(0);
        for( List<? extends Comparable> readGroupKey : readGroupTable.keySet() ) {

            if(NUM_READ_GROUPS_TO_PROCESS == -1 || ++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS) {
                String readGroup = readGroupKey.get(0).toString();
                RecalDatum readGroupDatum = readGroupTable.get(readGroupKey);
                System.out.print("Writing out data tables for read group: " + readGroup + "\twith " + readGroupDatum.getNumObservations() + " observations"  );
                System.out.println("\tand aggregate residual error = " + String.format("%.3f", readGroupDatum.empiricalQualDouble(0) - readGroupDatum.getEstimatedQReported()));

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

                    // Loop through the covariate table looking for keys with matching read groups
                    // BUGBUG: hopefully rewrite this to be more efficient
                    for( List<? extends Comparable> covariateKey : dataManager.getCollapsedTable(iii).keySet() ) {
                        if( covariateKey.get(0).toString().equals(readGroup) ) {
                            output.print( covariateKey.get(1).toString() + "\t" );                              // Covariate
                            RecalDatum thisDatum = dataManager.getCollapsedTable(iii).get(covariateKey);
                            output.print( String.format("%.3f", thisDatum.getEstimatedQReported()) + "\t" );    // Qreported
                            output.print( String.format("%.3f", thisDatum.empiricalQualDouble(0)) + "\t" );     // Qempirical
                            output.print( thisDatum.getNumMismatches() + "\t" );                                // nMismatches
                            output.println( thisDatum.getNumObservations() );                                   // nBases
                        }
                    }

                    // Close the PrintStream
                    output.close();
                }
            } else {
                break;
            }

        }
    }

    private static void callRScripts() {

        int numReadGroups = 0;
        
        // for each read group
        for( List<? extends Comparable> readGroupList : dataManager.getCollapsedTable(0).keySet() ) {

            if(NUM_READ_GROUPS_TO_PROCESS == -1 || ++numReadGroups <= NUM_READ_GROUPS_TO_PROCESS) {

                String readGroup = readGroupList.get(0).toString();
                System.out.println("Analyzing read group: " + readGroup);

                // for each covariate
                for( int iii = 1; iii < requestedCovariates.size(); iii++ ) {
                    Covariate cov = requestedCovariates.get(iii);
                    try {
                        if( iii == 1 ) {
                            // Analyze reported quality
                            Process p = Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_residualError_QualityScoreCovariate.R" + " " +
                                        OUTPUT_DIR + readGroup + "." + cov.getClass().getSimpleName()+ ".dat" + " " +
                                        IGNORE_QSCORES_LESS_THAN); // The third argument is the Q scores that should be turned pink in the plot because they were ignored
                        } else { // Analyze all other covariates
                            Process p = Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_residualError_OtherCovariate.R" + " " +
                                        OUTPUT_DIR + readGroup + "." + cov.getClass().getSimpleName()+ ".dat" + " " +
                                        cov.getClass().getSimpleName().split("Covariate")[0]); // The third argument which is the name of the covariate to make the plots look nice
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                        System.exit(-1);
                    }
                }
            } else {
                break;
            }
        }
    }
}

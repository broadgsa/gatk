/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Utility methods to facilitate on-the-fly base quality score recalibration.
 * 
 * User: rpoplin
 * Date: 2/4/12
 */

public class BaseRecalibration {

    private RecalDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>(); // List of covariates to be used in this calculation
    public static final Pattern COMMENT_PATTERN = Pattern.compile("^#.*");
    public static final Pattern COVARIATE_PATTERN = Pattern.compile("^ReadGroup,QualityScore,.*");
    public static final String EOF_MARKER = "EOF";
    private static final int MAX_QUALITY_SCORE = 65; //BUGBUG: what value to use here?
    private NestedHashMap qualityScoreByFullCovariateKey = new NestedHashMap(); // Caches the result of performSequentialQualityCalculation(...) for all sets of covariate values.

    public BaseRecalibration( final File RECAL_FILE ) {
        // Get a list of all available covariates
        final List<Class<? extends Covariate>> classes = new PluginManager<Covariate>(Covariate.class).getPlugins();

        int lineNumber = 0;
        boolean foundAllCovariates = false;

        // Read in the data from the csv file and populate the data map and covariates list
        boolean sawEOF = false;
        try {
            for ( String line : new XReadLines(RECAL_FILE) ) {
                lineNumber++;
                if ( EOF_MARKER.equals(line) ) {
                    sawEOF = true;
                } else if( COMMENT_PATTERN.matcher(line).matches() )  {
                    ; // Skip over the comment lines, (which start with '#')
                }
                // Read in the covariates that were used from the input file
                else if( COVARIATE_PATTERN.matcher(line).matches() ) { // The line string is either specifying a covariate or is giving csv data
                    if( foundAllCovariates ) {
                        throw new UserException.MalformedFile( RECAL_FILE, "Malformed input recalibration file. Found covariate names intermingled with data in file: " + RECAL_FILE );
                    } else { // Found the covariate list in input file, loop through all of them and instantiate them
                        String[] vals = line.split(",");
                        for( int iii = 0; iii < vals.length - 4; iii++ ) { // There are n-4 covariates. The last four items are ErrorModel, nObservations, nMismatch, and Qempirical
                            boolean foundClass = false;
                            for( Class<?> covClass : classes ) {
                                if( (vals[iii] + "Covariate").equalsIgnoreCase( covClass.getSimpleName() ) ) {
                                    foundClass = true;
                                    try {
                                        Covariate covariate = (Covariate)covClass.newInstance();
                                        requestedCovariates.add( covariate );
                                    } catch (Exception e) {
                                        throw new DynamicClassResolutionException(covClass, e);
                                    }

                                }
                            }

                            if( !foundClass ) {
                                throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration file. The requested covariate type (" + (vals[iii] + "Covariate") + ") isn't a valid covariate option." );
                            }
                        }
                    }

                } else { // Found a line of data
                    if( !foundAllCovariates ) {
                        foundAllCovariates = true;

                        // At this point all the covariates should have been found and initialized
                        if( requestedCovariates.size() < 2 ) {
                            throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Covariate names can't be found in file: " + RECAL_FILE );
                        }

                        final boolean createCollapsedTables = true;

                        // Initialize any covariate member variables using the shared argument collection
                        RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
                        for( Covariate cov : requestedCovariates ) {
                            cov.initialize( RAC );
                        }
                        // Initialize the data hashMaps
                        dataManager = new RecalDataManager( createCollapsedTables, requestedCovariates.size() );

                    }
                    addCSVData(RECAL_FILE, line); // Parse the line and add the data to the HashMap
                }
            }

        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(RECAL_FILE, "Can not find input file", e);
        } catch ( NumberFormatException e ) {
            throw new UserException.MalformedFile(RECAL_FILE, "Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
        }

        if ( !sawEOF ) {
            final String errorMessage = "No EOF marker was present in the recal covariates table; this could mean that the file is corrupted or was generated with an old version of the CountCovariates tool.";
            throw new UserException.MalformedFile(RECAL_FILE, errorMessage);
        }

        if( dataManager == null ) {
            throw new UserException.MalformedFile(RECAL_FILE, "Can't initialize the data manager. Perhaps the recal csv file contains no data?");
        }

        dataManager.generateEmpiricalQualities( 1, MAX_QUALITY_SCORE );
    }
    
    /**
     * For each covariate read in a value and parse it. Associate those values with the data itself (num observation and num mismatches)
     * @param line A line of CSV data read from the recalibration table data file
     */
    private void addCSVData(final File file, final String line) {
        final String[] vals = line.split(",");

        // Check if the data line is malformed, for example if the read group string contains a comma then it won't be parsed correctly
        if( vals.length != requestedCovariates.size() + 4 ) { // +4 because of ErrorModel, nObservations, nMismatch, and Qempirical
            throw new UserException.MalformedFile(file, "Malformed input recalibration file. Found data line with too many fields: " + line +
                    " --Perhaps the read group string contains a comma and isn't being parsed correctly.");
        }

        final Object[] key = new Object[requestedCovariates.size()];
        Covariate cov;
        int iii;
        for( iii = 0; iii < requestedCovariates.size(); iii++ ) {
            cov = requestedCovariates.get( iii );
            key[iii] = cov.getValue( vals[iii] );
        }
        final String modelString = vals[iii++];
        final RecalDataManager.BaseRecalibrationType errorModel = CovariateKeySet.getErrorModelFromString(modelString);

        // Create a new datum using the number of observations, number of mismatches, and reported quality score
        final RecalDatum datum = new RecalDatum( Long.parseLong( vals[iii] ), Long.parseLong( vals[iii + 1] ), Double.parseDouble( vals[1] ), 0.0 );
        // Add that datum to all the collapsed tables which will be used in the sequential calculation
        
        dataManager.addToAllTables( key, datum, QualityUtils.MIN_USABLE_Q_SCORE, errorModel ); //BUGBUG: used to be Q5 now is Q6, probably doesn't matter
    }
    
    public void recalibrateRead( final GATKSAMRecord read ) {

        //compute all covariate values for this read
        RecalDataManager.computeCovariates(read, requestedCovariates);
        final CovariateKeySet covariateKeySet = RecalDataManager.getAllCovariateValuesFor( read );

        for( final RecalDataManager.BaseRecalibrationType errorModel : RecalDataManager.BaseRecalibrationType.values() ) {
            final byte[] originalQuals = read.getBaseQualities( errorModel );
            final byte[] recalQuals = originalQuals.clone();

            // For each base in the read
            for( int offset = 0; offset < read.getReadLength(); offset++ ) {
        
                final Object[] fullCovariateKeyWithErrorMode = covariateKeySet.getKeySet(offset, errorModel);
                final Object[] fullCovariateKey = Arrays.copyOfRange(fullCovariateKeyWithErrorMode, 0, fullCovariateKeyWithErrorMode.length-1); // need to strip off the error mode which was appended to the list of covariates

                // BUGBUG: This caching seems to put the entire key set into memory which negates the benefits of storing the delta delta tables?
                //Byte qualityScore = (Byte) qualityScoreByFullCovariateKey.get(fullCovariateKeyWithErrorMode);
                //if( qualityScore == null ) {
                    final byte qualityScore = performSequentialQualityCalculation( errorModel, fullCovariateKey );
                //    qualityScoreByFullCovariateKey.put(qualityScore, fullCovariateKeyWithErrorMode);
                //}
        
                recalQuals[offset] = qualityScore;
            }
        
            preserveQScores( originalQuals, recalQuals ); // Overwrite the work done if original quality score is too low
            read.setBaseQualities( recalQuals, errorModel );
        }
    }

    /**
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     *   - calculate the global quality score shift across all data [DeltaQ]
     *   - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     *      -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     *   - The final shift equation is:
     *
     *      Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     * @param key The list of Comparables that were calculated from the covariates
     * @return A recalibrated quality score as a byte
     */
    private byte performSequentialQualityCalculation( final RecalDataManager.BaseRecalibrationType errorModel, final Object... key ) {

        final byte qualFromRead = (byte)Integer.parseInt(key[1].toString());
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] qualityScoreCollapsedKey = new Object[2];
        final Object[] covariateCollapsedKey = new Object[3];

        // The global quality shift (over the read group only)
        readGroupCollapsedKey[0] = key[0];
        final RecalDatum globalRecalDatum = ((RecalDatum)dataManager.getCollapsedTable(0, errorModel).get( readGroupCollapsedKey ));
        double globalDeltaQ = 0.0;
        if( globalRecalDatum != null ) {
            final double globalDeltaQEmpirical = globalRecalDatum.getEmpiricalQuality();
            final double aggregrateQReported = globalRecalDatum.getEstimatedQReported();
            globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
        }

        // The shift in quality between reported and empirical
        qualityScoreCollapsedKey[0] = key[0];
        qualityScoreCollapsedKey[1] = key[1];
        final RecalDatum qReportedRecalDatum = ((RecalDatum)dataManager.getCollapsedTable(1, errorModel).get( qualityScoreCollapsedKey ));
        double deltaQReported = 0.0;
        if( qReportedRecalDatum != null ) {
            final double deltaQReportedEmpirical = qReportedRecalDatum.getEmpiricalQuality();
            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        // The shift in quality due to each covariate by itself in turn
        double deltaQCovariates = 0.0;
        double deltaQCovariateEmpirical;
        covariateCollapsedKey[0] = key[0];
        covariateCollapsedKey[1] = key[1];
        for( int iii = 2; iii < key.length; iii++ ) {
            covariateCollapsedKey[2] =  key[iii]; // The given covariate
            final RecalDatum covariateRecalDatum = ((RecalDatum)dataManager.getCollapsedTable(iii, errorModel).get( covariateCollapsedKey ));
            if( covariateRecalDatum != null ) {
                deltaQCovariateEmpirical = covariateRecalDatum.getEmpiricalQuality();
                deltaQCovariates += ( deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported) );
            }
        }

        final double newQuality = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
        return QualityUtils.boundQual( (int)Math.round(newQuality), (byte)MAX_QUALITY_SCORE );
    }

    /**
     * Loop over the list of qualities and overwrite the newly recalibrated score to be the original score if it was less than some threshold
     * @param originalQuals The list of original base quality scores
     * @param recalQuals A list of the new recalibrated quality scores
     */
    private void preserveQScores( final byte[] originalQuals, final byte[] recalQuals ) {
        for( int iii = 0; iii < recalQuals.length; iii++ ) {
            if( originalQuals[iii] < QualityUtils.MIN_USABLE_Q_SCORE ) { //BUGBUG: used to be Q5 now is Q6, probably doesn't matter
                recalQuals[iii] = originalQuals[iii];
            }
        }
    }
}

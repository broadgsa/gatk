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
import org.broadinstitute.sting.utils.BitSetUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Utility methods to facilitate on-the-fly base quality score recalibration.
 *
 * User: rpoplin
 * Date: 2/4/12
 */

public class BaseRecalibration {

    private ArrayList<HashMap<BitSet, RecalDatum>> collapsedHashes = new ArrayList<HashMap<BitSet, RecalDatum>> (); // All the collapsed data tables

    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();                            // List of all covariates to be used in this calculation
    private final ArrayList<Covariate> requiredCovariates = new ArrayList<Covariate>();                             // List of required covariates to be used in this calculation
    private final ArrayList<Covariate> optionalCovariates = new ArrayList<Covariate>();                             // List of optional covariates to be used in this calculation
    
    public static final Pattern REQUIRED_COVARIATE_PATTERN = Pattern.compile("^# Required Covariates.*");
    public static final Pattern OPTIONAL_COVARIATE_PATTERN = Pattern.compile("^# Optional Covariates.*");
    public static final String EOF_MARKER = "EOF";

    private static final byte SMOOTHING_CONSTANT = 1;

    ArrayList<BQSRKeyManager> keyManagers = new ArrayList<BQSRKeyManager>();

    public BaseRecalibration(final File RECAL_FILE) {
        // Get a list of all available covariates
        final List<Class<? extends Covariate>> classes = new PluginManager<Covariate>(Covariate.class).getPlugins();
        RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection(); // todo -- initialize with the parameters from the csv file!

        int lineNumber = 0;

        boolean foundRequiredCovariates = false;
        boolean foundOptionalCovariates = false;
        boolean initializedKeyManagers = false;

        // Read in the data from the csv file and populate the data map and covariates list
        boolean sawEOF = false;
        try {
            for (String line : new XReadLines(RECAL_FILE)) {
                lineNumber++;

                sawEOF = EOF_MARKER.equals(line);
                if (sawEOF)
                    break;
                
                boolean requiredCovariatesLine = REQUIRED_COVARIATE_PATTERN.matcher(line).matches(); 
                boolean optionalCovariatesLine = OPTIONAL_COVARIATE_PATTERN.matcher(line).matches();
                
                if (requiredCovariatesLine && foundRequiredCovariates)
                    throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Duplicate required covariates line");

                if (optionalCovariatesLine && foundOptionalCovariates)
                    throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Duplicate optional covariates line");
                
                if (optionalCovariatesLine && !foundRequiredCovariates)
                    throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Optional covariates reported before Required covariates");
                
                if (requiredCovariatesLine || optionalCovariatesLine) {
                    String [] covariateNames = line.split(": ")[1].split(",");                                          // take the second half of the string (past the ":") and split it by "," to get the list of required covariates
                    
                    List<Covariate> covariateList = requiredCovariatesLine ? requiredCovariates : optionalCovariates;   // set the appropriate covariate list to update
                                
                    for (String covariateName : covariateNames) {
                        boolean foundClass = false;
                        for (Class<?> covClass : classes) {
                            if ((covariateName + "Covariate").equalsIgnoreCase(covClass.getSimpleName())) {
                                foundClass = true;
                                try {
                                    Covariate covariate = (Covariate) covClass.newInstance();
                                    covariate.initialize(RAC);
                                    requestedCovariates.add(covariate);
                                    covariateList.add(covariate);
                                } catch (Exception e) {
                                    throw new DynamicClassResolutionException(covClass, e);
                                }
                            }
                        }
                        if (!foundClass) 
                            throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration file. The requested covariate type (" + (covariateName + "Covariate") + ") isn't a valid covariate option.");
                    }
                    foundRequiredCovariates = foundRequiredCovariates || requiredCovariatesLine;
                    foundOptionalCovariates = foundOptionalCovariates || optionalCovariatesLine;
                }                

                else if (!line.startsWith("#")) {                                                                       // if this is not a comment line that we don't care about, it is DATA!
                    if (!foundRequiredCovariates || !foundOptionalCovariates)                                          // At this point all the covariates should have been found and initialized
                        throw new UserException.MalformedFile(RECAL_FILE, "Malformed input recalibration csv file. Covariate names can't be found in file: " + RECAL_FILE);

                    if (!initializedKeyManagers) {
                        ArrayList<Covariate> emptyList = new ArrayList<Covariate>(0);
                        ArrayList<Covariate> requiredCovariatesUpToThis = new ArrayList<Covariate>();                           // Initialize one key manager for each table of required covariate
                        for (Covariate covariate : requiredCovariates) {                                                // Every required covariate table includes all preceding required covariates (e.g. RG ; RG,Q )
                            requiredCovariatesUpToThis.add(covariate);
                            keyManagers.add(new BQSRKeyManager(requiredCovariatesUpToThis, emptyList));
                        }
                        keyManagers.add(new BQSRKeyManager(requiredCovariates, optionalCovariates));                   // One master key manager for the collapsed tables
                        
                        initializedKeyManagers = true;
                    }
                    addCSVData(RECAL_FILE, line);                                                                       // Parse the line and add the data to the HashMap
                }
            }

        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(RECAL_FILE, "Can not find input file", e);
        } catch (NumberFormatException e) {
            throw new UserException.MalformedFile(RECAL_FILE, "Error parsing recalibration data at line " + lineNumber + ". Perhaps your table was generated by an older version of CovariateCounterWalker.");
        }

        if (!sawEOF) {
            final String errorMessage = "No EOF marker was present in the recal covariates table; this could mean that the file is corrupted or was generated with an old version of the CountCovariates tool.";
            throw new UserException.MalformedFile(RECAL_FILE, errorMessage);
        }

        generateEmpiricalQualities(SMOOTHING_CONSTANT);
    }


    /**
     * For each covariate read in a value and parse it. Associate those values with the data itself (num observation and num mismatches)
     *
     * @param file The CSV file we read the line from (for exception throwing purposes)
     * @param line A line of CSV data read from the recalibration table data file
     */
    private void addCSVData(final File file, final String line) {
        final String[] vals = line.split(",");
        boolean hasOptionalCovariates = optionalCovariates.size() > 0;                                  // Do we have optional covariates in this key?
        int addOptionalCovariates = hasOptionalCovariates ? 2 : 0;                                      // If we have optional covariates at all, add two to the size of the array (to acommodate the covariate and the id)
        final Object[] key = new Object[requiredCovariates.size() + addOptionalCovariates + 1];         // Reserve enough space for the required covariates, optional covariate, id and eventType

        int indexCovariateValue = key.length - 3;                                                       // In the order of keys, the optional covariate comes right after the required covariates
        int indexCovariateID = key.length - 2;                                                          // followed by the covariate ID
        int indexEventType = key.length - 1;                                                            // and the event type

        addKeysToArray(key, vals, requiredCovariates, 0);                                               // Add the required covariates keys

        if (hasOptionalCovariates) {
            key[indexCovariateID] = Short.parseShort(vals[indexCovariateID]);                           // Add the optional covariate ID
            Covariate covariate = optionalCovariates.get((Short) key[indexCovariateID]);                // Get the covariate object for this ID
            key[indexCovariateValue] = covariate.getValue(vals[indexCovariateValue]);                   // Add the optional covariate value, given the ID
        }
        key[indexEventType] = EventType.eventFrom(vals[indexEventType]);                                // Add the event type 

        int datumIndex = key.length;                                                                    // The recal datum starts at the end of the key (after the event type)
        long count  = Long.parseLong(vals[datumIndex]);                                                 // Number of observations                                                 
        long errors = Long.parseLong(vals[datumIndex + 1]);                                             // Number of errors observed
        double reportedQual = Double.parseDouble(vals[1]);                                              // The reported Q score --> todo -- I don't like having the Q score hard coded in vals[1]. Generalize it!
        final RecalDatum datum = new RecalDatum(count, errors, reportedQual, 0.0);                      // Create a new datum using the number of observations, number of mismatches, and reported quality score

        addToAllTables(key, datum);                                                                     // Add that datum to all the collapsed tables which will be used in the sequential calculation
    }
    
    /**
     * Add the given mapping to all of the collapsed hash tables
     *
     * @param key       The list of comparables that is the key for this mapping
     * @param fullDatum The RecalDatum which is the data for this mapping
     */
    private void addToAllTables(final Object[] key, final RecalDatum fullDatum) {    
        int nHashes = requiredCovariates.size();                                                        // We will always need one hash per required covariate
        if (optionalCovariates.size() > 0)                                                              // If we do have optional covariates
            nHashes +=  1;                                                                              // we will need one extra hash table with the optional covariate encoded in the key set on top of the required covariates
        
        
        for (int hashIndex = 0; hashIndex < nHashes; hashIndex++) {
            HashMap<BitSet, RecalDatum> table;                                                          // object to hold the hash table we are going to manipulate
            if (hashIndex >= collapsedHashes.size()) {                                                  // if we haven't yet created the collapsed hash table for this index, create it now!
                table = new HashMap<BitSet, RecalDatum>();
                collapsedHashes.add(table);                                                             // Because this is the only place where we add tables to the ArrayList, they will always be in the order we want.
            }
            else
                table = collapsedHashes.get(hashIndex);                                                 // if the table has been previously created, just assign it to the "table" object for manipulation

            int copyTo = hashIndex + 1;                                                                 // this will copy the covariates up to the index of the one we are including now (1 for RG, 2 for QS,...)
            if (copyTo > requiredCovariates.size())                                                     // only in the case where we have optional covariates we need to increase the size of the array
                copyTo = requiredCovariates.size() + 2;                                                 // if we have optional covarites, add the optional covariate and it's id to the size of the key
            Object[] tableKey = new Object[copyTo + 1];                                                 // create a new array that will hold as many keys as hashIndex (1 for RG hash, 2 for QualityScore hash, 3 for covariate hash plus the event type
            System.arraycopy(key, 0, tableKey, 0, copyTo);                                              // copy the keys for the corresponding covariates into the tableKey.
            tableKey[tableKey.length-1] = key[key.length - 1];                                          // add the event type. The event type is always the last key, on both key sets.
            
            BitSet hashKey = keyManagers.get(hashIndex).bitSetFromKey(tableKey);                        // Add bitset key with fullDatum to the appropriate hash
            RecalDatum datum = table.get(hashKey);
            if (datum == null)
                datum = fullDatum;
            else if (hashIndex == 0)                                                                    // Special case for the ReadGroup covariate
                datum.combine(fullDatum);
            else
                datum.increment(fullDatum);
            table.put(hashKey, datum);
        }
    }


    /**
     * Loop over all the collapsed tables and turn the recalDatums found there into an empirical quality score
     * that will be used in the sequential calculation in TableRecalibrationWalker
     *
     * @param smoothing The smoothing parameter that goes into empirical quality score calculation
     */
    private void generateEmpiricalQualities(final int smoothing) {
        for (final HashMap<BitSet, RecalDatum> table : collapsedHashes)
            for (final RecalDatum datum : table.values())
                datum.calcCombinedEmpiricalQuality(smoothing, QualityUtils.MAX_QUAL_SCORE);
    }




    public void recalibrateRead(final GATKSAMRecord read) {           
        //compute all covariate values for this read
        RecalDataManager.computeCovariates(read, requestedCovariates);
        final ReadCovariates readCovariates = RecalDataManager.covariateKeySetFrom(read);

        for (final EventType errorModel : EventType.values()) {
            final byte[] originalQuals = read.getBaseQualities(errorModel);
            final byte[] recalQuals = originalQuals.clone();

            // For each base in the read
            for (int offset = 0; offset < read.getReadLength(); offset++) {
                final BitSet[] keySet = readCovariates.getKeySet(offset, errorModel);
                final byte qualityScore = performSequentialQualityCalculation(keySet, errorModel);
                recalQuals[offset] = qualityScore;
            }

            preserveQScores(originalQuals, recalQuals); // Overwrite the work done if original quality score is too low
            read.setBaseQualities(recalQuals, errorModel);
        }
    }


    
    /**
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     *
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     * 
     * todo -- I extremely dislike the way all this math is hardcoded... should rethink the data structures for this method in particular.
     *
     * @param key The list of Comparables that were calculated from the covariates
     * @param errorModel the event type
     * @return A recalibrated quality score as a byte
     */
    private byte performSequentialQualityCalculation(BitSet[] key, EventType errorModel) {
        final byte qualFromRead = (byte) BitSetUtils.shortFrom(key[1]);
               
        final int readGroupKeyIndex = 0;
        final int qualKeyIndex = 1;
        final int covariatesKeyIndex = 2;
                    
        // The global quality shift (over the read group only)
        List<BitSet> bitKeys = keyManagers.get(readGroupKeyIndex).bitSetsFromAllKeys(key, errorModel);
        if (bitKeys.size() > 1)
            throw new ReviewedStingException("There should only be one key for the RG collapsed table, something went wrong here");
        
        final RecalDatum globalRecalDatum = collapsedHashes.get(readGroupKeyIndex).get(bitKeys.get(0));
        double globalDeltaQ = 0.0;
        if (globalRecalDatum != null) {
            final double globalDeltaQEmpirical = globalRecalDatum.getEmpiricalQuality();
            final double aggregrateQReported = globalRecalDatum.getEstimatedQReported();
            globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
        }

        // The shift in quality between reported and empirical
        bitKeys = keyManagers.get(qualKeyIndex).bitSetsFromAllKeys(key, errorModel);
        if (bitKeys.size() > 1)
            throw new ReviewedStingException("There should only be one key for the Qual collapsed table, something went wrong here");

        final RecalDatum qReportedRecalDatum = collapsedHashes.get(qualKeyIndex).get(bitKeys.get(0));
        double deltaQReported = 0.0;
        if (qReportedRecalDatum != null) {
            final double deltaQReportedEmpirical = qReportedRecalDatum.getEmpiricalQuality();
            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        // The shift in quality due to each covariate by itself in turn
        bitKeys = keyManagers.get(covariatesKeyIndex).bitSetsFromAllKeys(key, errorModel);
        double deltaQCovariates = 0.0;
        double deltaQCovariateEmpirical;
        for (BitSet k : bitKeys) {
            final RecalDatum covariateRecalDatum = collapsedHashes.get(covariatesKeyIndex).get(k);
            if (covariateRecalDatum != null) {
                deltaQCovariateEmpirical = covariateRecalDatum.getEmpiricalQuality();
                deltaQCovariates += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
            }
        }

        final double newQuality = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
        return QualityUtils.boundQual((int) Math.round(newQuality), QualityUtils.MAX_QUAL_SCORE);
    }

    /**
     * Loop over the list of qualities and overwrite the newly recalibrated score to be the original score if it was less than some threshold
     *
     * @param originalQuals The list of original base quality scores
     * @param recalQuals    A list of the new recalibrated quality scores
     */
    private void preserveQScores(final byte[] originalQuals, final byte[] recalQuals) {
        for (int iii = 0; iii < recalQuals.length; iii++) {
            if (originalQuals[iii] < QualityUtils.MIN_USABLE_Q_SCORE) { //BUGBUG: used to be Q5 now is Q6, probably doesn't matter
                recalQuals[iii] = originalQuals[iii];
            }
        }
    }

    /**
     * Shared functionality to add keys
     * 
     * @param array         the target array we are creating the keys in
     * @param keys          the actual keys we're using as a source
     * @param covariateList the covariate list to loop through
     * @param keyIndex      the index in the keys and the arrays objects to run from
     */
    private void addKeysToArray(final Object[] array, final String[] keys, List<Covariate> covariateList, int keyIndex) {
        for (Covariate covariate : covariateList) {
            array[keyIndex] = covariate.getValue(keys[keyIndex]);        
            keyIndex++;
        }
    }
}

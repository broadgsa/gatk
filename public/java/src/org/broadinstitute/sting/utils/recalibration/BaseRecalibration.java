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

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.bqsr.*;
import org.broadinstitute.sting.gatk.walkers.recalibration.EmpiricalQual;
import org.broadinstitute.sting.utils.BitSetUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;
import java.util.*;

/**
 * Utility methods to facilitate on-the-fly base quality score recalibration.
 *
 * User: carneiro and rpoplin
 * Date: 2/4/12
 */

public class BaseRecalibration {
    private List<Byte> qualQuantizationMap;                                                                             // histogram containing the map for qual quantization (calculated after recalibration is done)
    private LinkedHashMap<BQSRKeyManager, Map<BitSet, EmpiricalQual>> keysAndTablesMap;                                 // quick access reference to the read group table and its key manager

    private ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();                                      // list of all covariates to be used in this calculation

    private static String UNRECOGNIZED_REPORT_TABLE_EXCEPTION = "Unrecognized table. Did you add an extra required covariate? This is a hard check that needs propagate through the code";
    private static String TOO_MANY_KEYS_EXCEPTION = "There should only be one key for the RG collapsed table, something went wrong here";

    /**
     * Should ALWAYS use the constructor with the GATK Report file 
     */
    private BaseRecalibration() {}

    /**
     * Constructor using a GATK Report file
     * 
     * @param RECAL_FILE a GATK Report file containing the recalibration information
     */
    public BaseRecalibration(final File RECAL_FILE) {
        GATKReport report = new GATKReport(RECAL_FILE);

        GATKReportTable argumentTable = report.getTable(RecalDataManager.ARGUMENT_REPORT_TABLE_TITLE);
        RecalibrationArgumentCollection RAC = initializeArgumentCollectionTable(argumentTable);

        GATKReportTable quantizedTable = report.getTable(RecalDataManager.QUANTIZED_REPORT_TABLE_TITLE);
        qualQuantizationMap = initializeQuantizationTable(quantizedTable);

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalDataManager.initializeCovariates(RAC);       // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getFirst();
        ArrayList<Covariate> optionalCovariates = covariates.getSecond();
        requestedCovariates.addAll(requiredCovariates);                                                                 // add all required covariates to the list of requested covariates
        requestedCovariates.addAll(optionalCovariates);                                                                 // add all optional covariates to the list of requested covariates

        for (Covariate cov : requestedCovariates) 
            cov.initialize(RAC);                                                                                        // initialize any covariate member variables using the shared argument collection
        
        keysAndTablesMap = new LinkedHashMap<BQSRKeyManager, Map<BitSet, EmpiricalQual>>();
        ArrayList<Covariate> requiredCovariatesToAdd = new ArrayList<Covariate>(requiredCovariates.size());                     // incrementally add the covariates to create the recal tables with 1, 2 and 3 covariates.
        ArrayList<Covariate> optionalCovariatesToAdd = new ArrayList<Covariate>();                                              // initialize an empty array of optional covariates to create the first few tables
        for (Covariate covariate : requiredCovariates) {
            requiredCovariatesToAdd.add(covariate);
            final Map<BitSet, EmpiricalQual> table;                                                                     // initializing a new recal table for each required covariate (cumulatively)
            final BQSRKeyManager keyManager = new BQSRKeyManager(requiredCovariatesToAdd, optionalCovariatesToAdd);     // initializing it's corresponding key manager
            
            int nRequiredCovariates = requiredCovariatesToAdd.size();                                                   // the number of required covariates defines which table we are looking at (RG, QUAL or ALL_COVARIATES)
            if (nRequiredCovariates == 1) {                                                                             // if there is only one required covariate, this is the read group table
                final GATKReportTable reportTable = report.getTable(RecalDataManager.READGROUP_REPORT_TABLE_TITLE);
                table = parseReadGroupTable(keyManager, reportTable);
            }
            else if (nRequiredCovariates == 2 && optionalCovariatesToAdd.isEmpty()) {                                   // when we have both required covariates and no optional covariates we're at the QUAL table
                final GATKReportTable reportTable = report.getTable(RecalDataManager.QUALITY_SCORE_REPORT_TABLE_TITLE);
                table = parseQualityScoreTable(keyManager, reportTable);
            }
            else
                throw new ReviewedStingException(UNRECOGNIZED_REPORT_TABLE_EXCEPTION);

            keysAndTablesMap.put(keyManager, table);                                                                    // adding the pair key+table to the map
        }


        final BQSRKeyManager keyManager = new BQSRKeyManager(requiredCovariates, optionalCovariates);                   // initializing it's corresponding key manager
        final GATKReportTable reportTable = report.getTable(RecalDataManager.ALL_COVARIATES_REPORT_TABLE_TITLE);
        final Map<BitSet, EmpiricalQual> table = parseAllCovariatesTable(keyManager, reportTable);
        keysAndTablesMap.put(keyManager, table);                                                                        // adding the pair table+key to the map
    }


    /**
     * Compiles the list of keys for the Covariates table and uses the shared parsing utility to produce the actual table
     *
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key. 
     */
    private Map<BitSet, EmpiricalQual> parseAllCovariatesTable(BQSRKeyManager keyManager, GATKReportTable reportTable) {
        ArrayList<String> columnNamesOrderedList = new ArrayList<String>(5);
        columnNamesOrderedList.add(RecalDataManager.READGROUP_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.QUALITY_SCORE_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.COVARIATE_VALUE_SCORE_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.COVARIATE_NAME_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.EVENT_TYPE_COLUMN_NAME);
        return genericRecalTableParsing(keyManager, reportTable, columnNamesOrderedList);
    }

    /**
     *
     * Compiles the list of keys for the QualityScore table and uses the shared parsing utility to produce the actual table
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key. 
     */
    private Map<BitSet, EmpiricalQual> parseQualityScoreTable(BQSRKeyManager keyManager, GATKReportTable reportTable) {
        ArrayList<String> columnNamesOrderedList = new ArrayList<String>(3);
        columnNamesOrderedList.add(RecalDataManager.READGROUP_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.QUALITY_SCORE_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.EVENT_TYPE_COLUMN_NAME);
        return genericRecalTableParsing(keyManager, reportTable, columnNamesOrderedList);
    }

    /**
     * Compiles the list of keys for the ReadGroup table and uses the shared parsing utility to produce the actual table
     *
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key. 
     */
    private Map<BitSet, EmpiricalQual> parseReadGroupTable(BQSRKeyManager keyManager, GATKReportTable reportTable) {
        ArrayList<String> columnNamesOrderedList = new ArrayList<String>(2);
        columnNamesOrderedList.add(RecalDataManager.READGROUP_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.EVENT_TYPE_COLUMN_NAME);
        return genericRecalTableParsing(keyManager, reportTable, columnNamesOrderedList);
    }

    /**
     * Shared parsing functionality for all tables.
     * 
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @param columnNamesOrderedList a list of columns to read from the report table and build as key for this particular table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key. 
     */
    private Map<BitSet, EmpiricalQual> genericRecalTableParsing(BQSRKeyManager keyManager, GATKReportTable reportTable, ArrayList<String> columnNamesOrderedList) {
        Map<BitSet, EmpiricalQual> result = new HashMap<BitSet, EmpiricalQual>(reportTable.getNumRows()*2);

        for (Object primaryKey : reportTable.getPrimaryKeys()) {
            int nKeys = columnNamesOrderedList.size();
            Object [] keySet = new Object[nKeys];
            for (int i = 0; i < nKeys; i++)
                keySet[i] = reportTable.get(primaryKey, columnNamesOrderedList.get(i));                                 // all these objects are okay in String format, the key manager will handle them correctly (except for the event type (see below)
            keySet[keySet.length-1] = EventType.eventFrom((String) keySet[keySet.length-1]);                            // the last key is always the event type. We convert the string ("M", "I" or "D") to an enum object (necessary for the key manager).
            BitSet bitKey = keyManager.bitSetFromKey(keySet);

            double estimatedQReported = (Double) reportTable.get(primaryKey, RecalDataManager.ESTIMATED_Q_REPORTED_COLUMN_NAME);
            double empiricalQuality = (Double) reportTable.get(primaryKey, RecalDataManager.EMPIRICAL_QUALITY_COLUMN_NAME);
            EmpiricalQual empiricalQual = new EmpiricalQual(estimatedQReported, empiricalQuality);

            result.put(bitKey, empiricalQual);
        }
        return result;
    }

    /**
     * Parses the quantization table from the GATK Report and turns it into a map of original => quantized quality scores
     * 
     * @param table the GATKReportTable containing the quantization mappings
     * @return an ArrayList with the quantization mappings from 0 to MAX_QUAL_SCORE
     */
    private List<Byte> initializeQuantizationTable(GATKReportTable table) {
        Byte[] result = new Byte[QualityUtils.MAX_QUAL_SCORE + 1];
        for (Object primaryKey : table.getPrimaryKeys()) {
            Object value = table.get(primaryKey, RecalDataManager.QUANTIZED_VALUE_COLUMN_NAME);
            byte originalQual = Byte.parseByte(primaryKey.toString());
            byte quantizedQual = Byte.parseByte(value.toString());
            result[originalQual] = quantizedQual;
        }
        return Arrays.asList(result);
    }

    /**
     * Parses the arguments table from the GATK Report and creates a RAC object with the proper initialization values
     *
     * @param table the GATKReportTable containing the arguments and its corresponding values
     * @return a RAC object properly initialized with all the objects in the table
     */
    private RecalibrationArgumentCollection initializeArgumentCollectionTable(GATKReportTable table) {
        RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        for (Object primaryKey : table.getPrimaryKeys()) {
            Object value = table.get(primaryKey, RecalDataManager.ARGUMENT_VALUE_COLUMN_NAME);
            if (value.equals("null"))
                value = null;                                                                                           // generic translation of null values that were printed out as strings | todo -- add this capability to the GATKReport
            
            if (primaryKey.equals("covariate") && value != null) 
                    RAC.COVARIATES = value.toString().split(",");            

            else if (primaryKey.equals("standard_covs"))
                RAC.USE_STANDARD_COVARIATES = Boolean.parseBoolean((String) value);
            
            else if (primaryKey.equals("solid_recal_mode"))
                RAC.SOLID_RECAL_MODE = RecalDataManager.SOLID_RECAL_MODE.recalModeFromString((String) value);

            else if (primaryKey.equals("solid_nocall_strategy"))
                RAC.SOLID_NOCALL_STRATEGY = RecalDataManager.SOLID_NOCALL_STRATEGY.nocallStrategyFromString((String) value);

            else if (primaryKey.equals("mismatches_context_size"))
                RAC.MISMATCHES_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (primaryKey.equals("insertions_context_size"))
                RAC.INSERTIONS_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (primaryKey.equals("deletions_context_size"))
                RAC.DELETIONS_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (primaryKey.equals("mismatches_default_quality"))
                RAC.MISMATCHES_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (primaryKey.equals("insertions_default_quality"))
                RAC.INSERTIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (primaryKey.equals("deletions_default_quality"))
                RAC.DELETIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (primaryKey.equals("low_quality_tail"))
                RAC.LOW_QUAL_TAIL = Byte.parseByte((String) value);

            else if (primaryKey.equals("default_platform"))
                RAC.DEFAULT_PLATFORM = (String) value;

            else if (primaryKey.equals("force_platform"))
                RAC.FORCE_PLATFORM = (String) value;

            else if (primaryKey.equals("quantizing_levels"))
                RAC.QUANTIZING_LEVELS = Integer.parseInt((String) value);
        }
        
        return RAC;
    }
       
    /**
     * Recalibrates the base qualities of a read
     *
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     *
     * @param read the read to recalibrate
     */
    public void recalibrateRead(final GATKSAMRecord read) {
        final ReadCovariates readCovariates = RecalDataManager.computeCovariates(read, requestedCovariates);    // compute all covariates for the read
        for (final EventType errorModel : EventType.values()) {                                                 // recalibrate all three quality strings
            final byte[] originalQuals = read.getBaseQualities(errorModel);
            final byte[] recalQuals = originalQuals.clone();

            for (int offset = 0; offset < read.getReadLength(); offset++) {                                     // recalibrate all bases in the read
                byte qualityScore = originalQuals[offset];

                if (qualityScore > QualityUtils.MIN_USABLE_Q_SCORE) {                                           // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                    final BitSet[] keySet = readCovariates.getKeySet(offset, errorModel);                       // get the keyset for this base using the error model
                    qualityScore = performSequentialQualityCalculation(keySet, errorModel);                     // recalibrate the base
                }
                recalQuals[offset] = qualityScore;
            }
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
     * @param key        The list of Comparables that were calculated from the covariates
     * @param errorModel the event type
     * @return A recalibrated quality score as a byte
     */
    private byte performSequentialQualityCalculation(BitSet[] key, EventType errorModel) {
        final byte qualFromRead = (byte) BitSetUtils.shortFrom(key[1]);

        double globalDeltaQ = 0.0;
        double deltaQReported = 0.0;
        double deltaQCovariates = 0.0;

        for (Map.Entry<BQSRKeyManager, Map<BitSet, EmpiricalQual>> mapEntry : keysAndTablesMap.entrySet()) {
            BQSRKeyManager keyManager = mapEntry.getKey();
            Map<BitSet, EmpiricalQual> table = mapEntry.getValue();

            switch(keyManager.getRequiredCovariates().size()) {
                case 1:                                                                                                 // this is the ReadGroup table                    
                    List<BitSet> bitKeys = keyManager.bitSetsFromAllKeys(key, errorModel);                              // calculate the shift in quality due to the read group
                    if (bitKeys.size() > 1)
                        throw new ReviewedStingException(TOO_MANY_KEYS_EXCEPTION);

                    final EmpiricalQual empiricalQualRG = table.get(bitKeys.get(0));
                    if (empiricalQualRG != null) {
                        final double globalDeltaQEmpirical = empiricalQualRG.getEmpiricalQuality();
                        final double aggregrateQReported = empiricalQualRG.getEstimatedQReported();
                        globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
                    }
                    break;
                case 2:
                    if (keyManager.getOptionalCovariates().isEmpty()) {                                                 // this is the QualityScore table
                        bitKeys = keyManager.bitSetsFromAllKeys(key, errorModel);                                       // calculate the shift in quality due to the reported quality score
                        if (bitKeys.size() > 1)
                            throw new ReviewedStingException(TOO_MANY_KEYS_EXCEPTION);

                        final EmpiricalQual empiricalQualQS = table.get(bitKeys.get(0));
                        if (empiricalQualQS != null) {
                            final double deltaQReportedEmpirical = empiricalQualQS.getEmpiricalQuality();
                            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
                        }
                    }
                    else {                                                                                              // this is the table with all the covariates                        
                        bitKeys = keyManager.bitSetsFromAllKeys(key, errorModel);                                       // calculate the shift in quality due to each covariate by itself in turn
                        for (BitSet k : bitKeys) {
                            final EmpiricalQual empiricalQualCO = table.get(k);
                            if (empiricalQualCO != null) {
                                double deltaQCovariateEmpirical = empiricalQualCO.getEmpiricalQuality();
                                deltaQCovariates += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
                            }
                        }
                    }
                    break;
                default:
                    throw new ReviewedStingException(UNRECOGNIZED_REPORT_TABLE_EXCEPTION);
            }
        }

        double recalibratedQual = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;                      // calculate the recalibrated qual using the BQSR formula 
        recalibratedQual = QualityUtils.boundQual((int) Math.round(recalibratedQual), QualityUtils.MAX_QUAL_SCORE);     // recalibrated quality is bound between 1 and MAX_QUAL

        return qualQuantizationMap.get((int) recalibratedQual);                                                         // return the quantized version of the recalibrated quality        
    }

}

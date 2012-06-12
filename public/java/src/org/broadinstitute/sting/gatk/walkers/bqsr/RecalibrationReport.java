package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * This class has all the static functionality for reading a recalibration report file into memory. 
 *
 * @author carneiro
 * @since 3/26/12
 */
public class RecalibrationReport {
    private QuantizationInfo quantizationInfo;                                                                          // histogram containing the counts for qual quantization (calculated after recalibration is done)
    private final LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> keysAndTablesMap;                                // quick access reference to the read group table and its key manager
    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();                                        // list of all covariates to be used in this calculation

    private final GATKReportTable argumentTable;                                                                              // keep the argument table untouched just for output purposes
    private final RecalibrationArgumentCollection RAC;                                                                        // necessary for quantizing qualities with the same parameter

    public RecalibrationReport(final File RECAL_FILE) {
        GATKReport report = new GATKReport(RECAL_FILE);

        argumentTable = report.getTable(RecalDataManager.ARGUMENT_REPORT_TABLE_TITLE);
        RAC = initializeArgumentCollectionTable(argumentTable);

        GATKReportTable quantizedTable = report.getTable(RecalDataManager.QUANTIZED_REPORT_TABLE_TITLE);
        quantizationInfo = initializeQuantizationTable(quantizedTable);

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalDataManager.initializeCovariates(RAC);       // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getFirst();
        ArrayList<Covariate> optionalCovariates = covariates.getSecond();
        requestedCovariates.addAll(requiredCovariates);                                                                 // add all required covariates to the list of requested covariates
        requestedCovariates.addAll(optionalCovariates);                                                                 // add all optional covariates to the list of requested covariates

        for (Covariate cov : requestedCovariates)
            cov.initialize(RAC);                                                                                        // initialize any covariate member variables using the shared argument collection

        keysAndTablesMap = new LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>>();
        ArrayList<Covariate> requiredCovariatesToAdd = new ArrayList<Covariate>(requiredCovariates.size());                     // incrementally add the covariates to create the recal tables with 1, 2 and 3 covariates.
        ArrayList<Covariate> optionalCovariatesToAdd = new ArrayList<Covariate>();                                              // initialize an empty array of optional covariates to create the first few tables
        for (Covariate covariate : requiredCovariates) {
            requiredCovariatesToAdd.add(covariate);
            final Map<Long, RecalDatum> table;                                                                          // initializing a new recal table for each required covariate (cumulatively)
            final BQSRKeyManager keyManager = new BQSRKeyManager(requiredCovariatesToAdd, optionalCovariatesToAdd);     // initializing it's corresponding key manager

            final int nRequiredCovariates = requiredCovariatesToAdd.size();                                             // the number of required covariates defines which table we are looking at (RG, QUAL or ALL_COVARIATES)
            final String UNRECOGNIZED_REPORT_TABLE_EXCEPTION = "Unrecognized table. Did you add an extra required covariate? This is a hard check.";
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
        final Map<Long, RecalDatum> table = parseAllCovariatesTable(keyManager, reportTable);
        keysAndTablesMap.put(keyManager, table);
    }

    protected RecalibrationReport(QuantizationInfo quantizationInfo, LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> keysAndTablesMap, GATKReportTable argumentTable, RecalibrationArgumentCollection RAC) {
        this.quantizationInfo = quantizationInfo;
        this.keysAndTablesMap = keysAndTablesMap;
        this.argumentTable = argumentTable;
        this.RAC = RAC;
    }

    /**
    * Combines two recalibration reports by adding all observations and errors
    *
    * Note: This method DOES NOT recalculate the empirical qualities and quantized qualities. You have to recalculate
    * them after combining. The reason for not calculating it is because this function is inteded for combining a
    * series of recalibration reports, and it only makes sense to calculate the empirical qualities and quantized
    * qualities after all the recalibration reports have been combined. Having the user recalculate when appropriate,
    * makes this method faster
    *
    * Note2: The empirical quality reported, however, is recalculated given its simplicity.
    *
    * @param other the recalibration report to combine with this one
    */
    public void combine(RecalibrationReport other) {
        Iterator<Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>>> thisIterator = keysAndTablesMap.entrySet().iterator();

        for (Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> otherEntry : other.getKeysAndTablesMap().entrySet()) {
            Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> thisEntry = thisIterator.next();

            final Map<Long, RecalDatum> thisTable = thisEntry.getValue();
            final BQSRKeyManager thisKeyManager = thisEntry.getKey();
            final BQSRKeyManager otherKeyManager = otherEntry.getKey();

            for (Map.Entry<Long, RecalDatum> otherTableEntry : otherEntry.getValue().entrySet()) {
                final RecalDatum otherDatum = otherTableEntry.getValue();
                final Long otherBitKey = otherTableEntry.getKey();
                final List<Object> otherObjectKey = otherKeyManager.keySetFrom(otherBitKey);
                
                final Long thisKey = thisKeyManager.longFromKey(otherObjectKey.toArray());
                final RecalDatum thisDatum = thisTable.get(thisKey);
                
                if (thisDatum == null)
                    thisTable.put(thisKey, otherDatum);
                else
                    thisDatum.combine(otherDatum);
            }            
        }
    }

    public QuantizationInfo getQuantizationInfo() {
        return quantizationInfo;
    }

    public LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> getKeysAndTablesMap() {
        return keysAndTablesMap;
    }

    public ArrayList<Covariate> getRequestedCovariates() {
        return requestedCovariates;
    }

    /**
     * Compiles the list of keys for the Covariates table and uses the shared parsing utility to produce the actual table
     *
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key.
     */
    private Map<Long, RecalDatum> parseAllCovariatesTable(BQSRKeyManager keyManager, GATKReportTable reportTable) {
        ArrayList<String> columnNamesOrderedList = new ArrayList<String>(5);
        columnNamesOrderedList.add(RecalDataManager.READGROUP_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.QUALITY_SCORE_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.COVARIATE_VALUE_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.COVARIATE_NAME_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.EVENT_TYPE_COLUMN_NAME);
        return genericRecalTableParsing(keyManager, reportTable, columnNamesOrderedList, false);
    }

    /**
     *
     * Compiles the list of keys for the QualityScore table and uses the shared parsing utility to produce the actual table
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key.
     */
    private Map<Long, RecalDatum> parseQualityScoreTable(BQSRKeyManager keyManager, GATKReportTable reportTable) {
        ArrayList<String> columnNamesOrderedList = new ArrayList<String>(3);
        columnNamesOrderedList.add(RecalDataManager.READGROUP_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.QUALITY_SCORE_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.EVENT_TYPE_COLUMN_NAME);
        return genericRecalTableParsing(keyManager, reportTable, columnNamesOrderedList, false);
    }

    /**
     * Compiles the list of keys for the ReadGroup table and uses the shared parsing utility to produce the actual table
     *
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key.
     */
    private Map<Long, RecalDatum> parseReadGroupTable(BQSRKeyManager keyManager, GATKReportTable reportTable) {
        ArrayList<String> columnNamesOrderedList = new ArrayList<String>(2);
        columnNamesOrderedList.add(RecalDataManager.READGROUP_COLUMN_NAME);
        columnNamesOrderedList.add(RecalDataManager.EVENT_TYPE_COLUMN_NAME);
        return genericRecalTableParsing(keyManager, reportTable, columnNamesOrderedList, true);
    }

    /**
     * Shared parsing functionality for all tables.
     *
     * @param keyManager             the key manager for this table
     * @param reportTable            the GATKReport table containing data for this table
     * @param columnNamesOrderedList a list of columns to read from the report table and build as key for this particular table
     * @return a lookup table indexed by bitsets containing the empirical quality and estimated quality reported for every key.
     */
    private Map<Long, RecalDatum> genericRecalTableParsing(BQSRKeyManager keyManager, GATKReportTable reportTable, ArrayList<String> columnNamesOrderedList, boolean hasEstimatedQReportedColumn) {
        final Map<Long, RecalDatum> result = new HashMap<Long, RecalDatum>(reportTable.getNumRows()*2);

        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final int nKeys = columnNamesOrderedList.size();
            final Object [] keySet = new Object[nKeys];
            for (int j = 0; j < nKeys; j++)
                keySet[j] = reportTable.get(i, columnNamesOrderedList.get(j));                                          // all these objects are okay in String format, the key manager will handle them correctly (except for the event type (see below)
            keySet[keySet.length-1] = EventType.eventFrom((String) keySet[keySet.length-1]);                            // the last key is always the event type. We convert the string ("M", "I" or "D") to an enum object (necessary for the key manager).
            final Long bitKey = keyManager.longFromKey(keySet);

            final long nObservations = (Long) reportTable.get(i, RecalDataManager.NUMBER_OBSERVATIONS_COLUMN_NAME);
            final long nErrors = (Long) reportTable.get(i, RecalDataManager.NUMBER_ERRORS_COLUMN_NAME);
            final double empiricalQuality = (Double) reportTable.get(i, RecalDataManager.EMPIRICAL_QUALITY_COLUMN_NAME);

            final double estimatedQReported = hasEstimatedQReportedColumn ?                                             // the estimatedQreported column only exists in the ReadGroup table
                (Double) reportTable.get(i, RecalDataManager.ESTIMATED_Q_REPORTED_COLUMN_NAME) :                        // we get it if we are in the read group table
                Byte.parseByte((String) reportTable.get(i, RecalDataManager.QUALITY_SCORE_COLUMN_NAME));                // or we use the reported quality if we are in any other table

            final RecalDatum recalDatum = new RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);
            result.put(bitKey, recalDatum);
        }
        return result;
    }

    /**
     * Parses the quantization table from the GATK Report and turns it into a map of original => quantized quality scores
     *
     * @param table the GATKReportTable containing the quantization mappings
     * @return an ArrayList with the quantization mappings from 0 to MAX_QUAL_SCORE
     */
    private QuantizationInfo initializeQuantizationTable(GATKReportTable table) {
        final Byte[] quals  = new Byte[QualityUtils.MAX_QUAL_SCORE + 1];
        final Long[] counts = new Long[QualityUtils.MAX_QUAL_SCORE + 1];
        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final byte originalQual = (byte)i;
            final Object quantizedObject = table.get(i, RecalDataManager.QUANTIZED_VALUE_COLUMN_NAME);
            final Object countObject = table.get(i, RecalDataManager.QUANTIZED_COUNT_COLUMN_NAME);
            final byte quantizedQual = Byte.parseByte(quantizedObject.toString());
            final long quantizedCount = Long.parseLong(countObject.toString());
            quals[originalQual] = quantizedQual;
            counts[originalQual] = quantizedCount;
        }
        return new QuantizationInfo(Arrays.asList(quals), Arrays.asList(counts));
    }

    /**
     * Parses the arguments table from the GATK Report and creates a RAC object with the proper initialization values
     *
     * @param table the GATKReportTable containing the arguments and its corresponding values
     * @return a RAC object properly initialized with all the objects in the table
     */
    private RecalibrationArgumentCollection initializeArgumentCollectionTable(GATKReportTable table) {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final String argument = table.get(i, "Argument").toString();
            Object value = table.get(i, RecalDataManager.ARGUMENT_VALUE_COLUMN_NAME);
            if (value.equals("null"))
                value = null;                                                                                           // generic translation of null values that were printed out as strings | todo -- add this capability to the GATKReport

            if (argument.equals("covariate") && value != null)
                RAC.COVARIATES = value.toString().split(",");

            else if (argument.equals("standard_covs"))
                RAC.USE_STANDARD_COVARIATES = Boolean.parseBoolean((String) value);

            else if (argument.equals("solid_recal_mode"))
                RAC.SOLID_RECAL_MODE = RecalDataManager.SOLID_RECAL_MODE.recalModeFromString((String) value);

            else if (argument.equals("solid_nocall_strategy"))
                RAC.SOLID_NOCALL_STRATEGY = RecalDataManager.SOLID_NOCALL_STRATEGY.nocallStrategyFromString((String) value);

            else if (argument.equals("mismatches_context_size"))
                RAC.MISMATCHES_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("insertions_context_size"))
                RAC.INSERTIONS_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("deletions_context_size"))
                RAC.DELETIONS_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("mismatches_default_quality"))
                RAC.MISMATCHES_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("insertions_default_quality"))
                RAC.INSERTIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("deletions_default_quality"))
                RAC.DELETIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("low_quality_tail"))
                RAC.LOW_QUAL_TAIL = Byte.parseByte((String) value);

            else if (argument.equals("default_platform"))
                RAC.DEFAULT_PLATFORM = (String) value;

            else if (argument.equals("force_platform"))
                RAC.FORCE_PLATFORM = (String) value;

            else if (argument.equals("quantizing_levels"))
                RAC.QUANTIZING_LEVELS = Integer.parseInt((String) value);

            else if (argument.equals("keep_intermediate_files"))
                RAC.KEEP_INTERMEDIATE_FILES = Boolean.parseBoolean((String) value);

            else if (argument.equals("no_plots"))
                RAC.NO_PLOTS = Boolean.parseBoolean((String) value);

            else if (argument.equals("recalibration_report"))
                RAC.recalibrationReport = (value == null) ? null : new File((String) value);
        }

        return RAC;
    }

    /**
     * this functionality avoids recalculating the empirical qualities, estimated reported quality
     * and quantization of the quality scores during every call of combine(). Very useful for the BQSRGatherer.
     */
    public void calculateEmpiricalAndQuantizedQualities() {
        for (Map<Long, RecalDatum> table : keysAndTablesMap.values())
            for (RecalDatum datum : table.values())
                datum.calcCombinedEmpiricalQuality();

        quantizationInfo = new QuantizationInfo(keysAndTablesMap, RAC.QUANTIZING_LEVELS);
    }

    public void output(PrintStream output) {
        RecalDataManager.outputRecalibrationReport(argumentTable, quantizationInfo, keysAndTablesMap, output);
    }

    public RecalibrationArgumentCollection getRAC() {
        return RAC;
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof RecalibrationReport))
            return false;
        RecalibrationReport other = (RecalibrationReport) o;
        if (this == o)
            return true;
        return isEqualTable(this.keysAndTablesMap, other.keysAndTablesMap);
    }

    private boolean isEqualTable(LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> t1, LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> t2) {
        if (t1.size() != t2.size())
            return false;

        final Iterator<Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>>> t1Iterator = t1.entrySet().iterator();
        final Iterator<Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>>> t2Iterator = t2.entrySet().iterator();

        while (t1Iterator.hasNext() && t2Iterator.hasNext()) {
            Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> t1MapEntry = t1Iterator.next();
            Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> t2MapEntry = t2Iterator.next();

            if (!(t1MapEntry.getKey().equals(t2MapEntry.getKey())))
                return false;

            final Map<Long, RecalDatum> table2 = t2MapEntry.getValue();
            for (Map.Entry<Long, RecalDatum> t1TableEntry : t1MapEntry.getValue().entrySet()) {
                final Long t1Key = t1TableEntry.getKey();
                if (!table2.containsKey(t1Key))
                    return false;
                final RecalDatum t1Datum = t1TableEntry.getValue();
                if (!t1Datum.equals(table2.get(t1Key)))
                    return false;
            }
        }
        return true;
    }
}

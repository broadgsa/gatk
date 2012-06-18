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
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
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
    private final static int MAXIMUM_RECALIBRATED_READ_LENGTH = 5000;
    private final ReadCovariates readCovariates;

    private final QuantizationInfo quantizationInfo;                                                                    // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final KeysAndTables keysAndTables;
    private final Covariate[] requestedCovariates;                                                                      // list of all covariates to be used in this calculation

    static class KeysAndTables {

        public enum Type {
            READ_GROUP_TABLE(0),
            QUALITY_SCORE_TABLE(1),
            OPTIONAL_COVARIATE_TABLE(2);

            private final int index;

            private Type(int index) {
                this.index = index;
            }
        }

        public final BQSRKeyManager[] managers = new BQSRKeyManager[Type.values().length];
        public final Map<Long, RecalDatum>[] tables = new Map[Type.values().length];

        public KeysAndTables(final Map<BQSRKeyManager, Map<Long, RecalDatum>> keysAndTablesMap) {
            for (Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> mapEntry : keysAndTablesMap.entrySet()) {
                Type type;
                if (mapEntry.getKey().getNumRequiredCovariates() == 1)
                    type = Type.READ_GROUP_TABLE;
                else if (mapEntry.getKey().getNumOptionalCovariates() == 0)
                    type = Type.QUALITY_SCORE_TABLE;
                else
                    type = Type.OPTIONAL_COVARIATE_TABLE;
                managers[type.index] = mapEntry.getKey();
                tables[type.index] = mapEntry.getValue();
            }
        }
    }

    /**
     * Constructor using a GATK Report file
     * 
     * @param RECAL_FILE a GATK Report file containing the recalibration information
     * @param quantizationLevels number of bins to quantize the quality scores
     */
    public BaseRecalibration(final File RECAL_FILE, int quantizationLevels) {
        RecalibrationReport recalibrationReport = new RecalibrationReport(RECAL_FILE);

        keysAndTables = new KeysAndTables(recalibrationReport.getKeysAndTablesMap());
        requestedCovariates = recalibrationReport.getRequestedCovariates();
        quantizationInfo = recalibrationReport.getQuantizationInfo();
        if (quantizationLevels == 0)                                                                                    // quantizationLevels == 0 means no quantization, preserve the quality scores
            quantizationInfo.noQuantization();
        else if (quantizationLevels > 0 && quantizationLevels != quantizationInfo.getQuantizationLevels())              // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wnats to use what's in the report.
            quantizationInfo.quantizeQualityScores(quantizationLevels);

        readCovariates = new ReadCovariates(MAXIMUM_RECALIBRATED_READ_LENGTH, requestedCovariates.length);
    }

    /**
     * This constructor only exists for testing purposes.
     *
     * @param quantizationInfo the quantization info object
     * @param keysAndTablesMap the map of key managers and recalibration tables
     * @param requestedCovariates the list of requested covariates
     */
    protected BaseRecalibration(final QuantizationInfo quantizationInfo, final LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> keysAndTablesMap, final Covariate[] requestedCovariates) {
        this.quantizationInfo = quantizationInfo;
        keysAndTables = new KeysAndTables(keysAndTablesMap);
        this.requestedCovariates = requestedCovariates;
        readCovariates = new ReadCovariates(MAXIMUM_RECALIBRATED_READ_LENGTH, requestedCovariates.length);
    }

    /**
     * Recalibrates the base qualities of a read
     *
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     *
     * @param read the read to recalibrate
     */
    public void recalibrateRead(final GATKSAMRecord read) {
        RecalDataManager.computeCovariates(read, requestedCovariates, readCovariates);                                  // compute all covariates for the read
        for (final EventType errorModel : EventType.values()) {                                                         // recalibrate all three quality strings
            final byte[] quals = read.getBaseQualities(errorModel);

            for (int offset = 0; offset < read.getReadLength(); offset++) {                                             // recalibrate all bases in the read
                final byte originalQualityScore = quals[offset];

                if (originalQualityScore >= QualityUtils.MIN_USABLE_Q_SCORE) {                                          // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                    final long[] keySet = readCovariates.getKeySet(offset, errorModel);                                 // get the keyset for this base using the error model
                    final byte recalibratedQualityScore = performSequentialQualityCalculation(keySet, errorModel);      // recalibrate the base
                    quals[offset] = recalibratedQualityScore;
                }
            }
            read.setBaseQualities(quals, errorModel);
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
    protected byte performSequentialQualityCalculation(final long[] key, final EventType errorModel) {

        final double globalDeltaQ = calculateGlobalDeltaQ(keysAndTables.managers[KeysAndTables.Type.READ_GROUP_TABLE.index], keysAndTables.tables[KeysAndTables.Type.READ_GROUP_TABLE.index], key, errorModel);
        final double deltaQReported = calculateDeltaQReported(keysAndTables.managers[KeysAndTables.Type.QUALITY_SCORE_TABLE.index], keysAndTables.tables[KeysAndTables.Type.QUALITY_SCORE_TABLE.index], key, errorModel, globalDeltaQ);
        final double deltaQCovariates = calculateDeltaQCovariates(keysAndTables.managers[KeysAndTables.Type.OPTIONAL_COVARIATE_TABLE.index], keysAndTables.tables[KeysAndTables.Type.OPTIONAL_COVARIATE_TABLE.index], key, errorModel, globalDeltaQ, deltaQReported);

        final byte qualFromRead = (byte)key[1];
        double recalibratedQual = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;                      // calculate the recalibrated qual using the BQSR formula
        recalibratedQual = QualityUtils.boundQual(MathUtils.fastRound(recalibratedQual), QualityUtils.MAX_RECALIBRATED_Q_SCORE);     // recalibrated quality is bound between 1 and MAX_QUAL

        return quantizationInfo.getQuantizedQuals().get((int) recalibratedQual);                                        // return the quantized version of the recalibrated quality
    }

    private double calculateGlobalDeltaQ(final BQSRKeyManager keyManager, final Map<Long, RecalDatum> table, final long[] key, final EventType errorModel) {
        double result = 0.0;

        final long masterKey = keyManager.createMasterKey(key, errorModel, -1);
        final RecalDatum empiricalQualRG = table.get(masterKey);
        if (empiricalQualRG != null) {
            final double globalDeltaQEmpirical = empiricalQualRG.getEmpiricalQuality();
            final double aggregrateQReported = empiricalQualRG.getEstimatedQReported();
            result = globalDeltaQEmpirical - aggregrateQReported;
        }

        return result;
    }

    private double calculateDeltaQReported(final BQSRKeyManager keyManager, final Map<Long, RecalDatum> table, final long[] key, final EventType errorModel, final double globalDeltaQ) {
        double result = 0.0;

        final long masterKey = keyManager.createMasterKey(key, errorModel, -1);
        final RecalDatum empiricalQualQS = table.get(masterKey);
        if (empiricalQualQS != null) {
            final double deltaQReportedEmpirical = empiricalQualQS.getEmpiricalQuality();
            final byte qualFromRead = (byte)key[1];
            result = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        return result;
    }

    private double calculateDeltaQCovariates(final BQSRKeyManager keyManager, final Map<Long, RecalDatum> table, final long[] key, final EventType errorModel, final double globalDeltaQ, final double deltaQReported) {
        double result = 0.0;

        final int numOptionalCovariates = keyManager.getNumOptionalCovariates();
        for (int i = 0; i < numOptionalCovariates; i++) {
            final long masterKey = keyManager.createMasterKey(key, errorModel, i);
            if (masterKey < 0)
                continue;
            final RecalDatum empiricalQualCO = table.get(masterKey);
            if (empiricalQualCO != null) {
                final double deltaQCovariateEmpirical = empiricalQualCO.getEmpiricalQuality();
                final byte qualFromRead = (byte)key[1];
                result += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
            }
        }
        return result;
    }
}

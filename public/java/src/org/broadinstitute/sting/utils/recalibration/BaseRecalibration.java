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

import net.sf.samtools.SAMTag;
import net.sf.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;

/**
 * Utility methods to facilitate on-the-fly base quality score recalibration.
 *
 * User: carneiro and rpoplin
 * Date: 2/4/12
 */

public class BaseRecalibration {
    private static Logger logger = Logger.getLogger(BaseRecalibration.class);
    private final static boolean TEST_CACHING = false;

    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final Covariate[] requestedCovariates; // list of all covariates to be used in this calculation

    private final boolean disableIndelQuals;
    private final int preserveQLessThan;
    private final boolean emitOriginalQuals;

    private final NestedIntegerArray<Double> globalDeltaQs;
    private final NestedIntegerArray<Double> deltaQReporteds;


    /**
     * Constructor using a GATK Report file
     * 
     * @param RECAL_FILE         a GATK Report file containing the recalibration information
     * @param quantizationLevels number of bins to quantize the quality scores
     * @param disableIndelQuals  if true, do not emit base indel qualities
     * @param preserveQLessThan  preserve quality scores less than this value
     */
    public BaseRecalibration(final File RECAL_FILE, final int quantizationLevels, final boolean disableIndelQuals, final int preserveQLessThan, final boolean emitOriginalQuals) {
        RecalibrationReport recalibrationReport = new RecalibrationReport(RECAL_FILE);

        recalibrationTables = recalibrationReport.getRecalibrationTables();
        requestedCovariates = recalibrationReport.getRequestedCovariates();
        quantizationInfo = recalibrationReport.getQuantizationInfo();
        if (quantizationLevels == 0) // quantizationLevels == 0 means no quantization, preserve the quality scores
            quantizationInfo.noQuantization();
        else if (quantizationLevels > 0 && quantizationLevels != quantizationInfo.getQuantizationLevels()) // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
            quantizationInfo.quantizeQualityScores(quantizationLevels);

        this.disableIndelQuals = disableIndelQuals;
        this.preserveQLessThan = preserveQLessThan;
        this.emitOriginalQuals = emitOriginalQuals;

        logger.info("Calculating cached tables...");

        //
        // Create a NestedIntegerArray<Double> that maps from rgKey x errorModel -> double,
        // where the double is the result of this calculation.  The entire calculation can
        // be done upfront, on initialization of this BaseRecalibration structure
        //
        final NestedIntegerArray<RecalDatum> byReadGroupTable = recalibrationTables.getReadGroupTable();
        globalDeltaQs = new NestedIntegerArray<Double>( byReadGroupTable.getDimensions() );
        logger.info("Calculating global delta Q table...");
        for ( NestedIntegerArray.Leaf<RecalDatum> leaf : byReadGroupTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[1];
            final double globalDeltaQ = calculateGlobalDeltaQ(rgKey, EventType.eventFrom(eventIndex));
            globalDeltaQs.put(globalDeltaQ, rgKey, eventIndex);
        }


        // The calculation of the deltaQ report is constant.  key[0] and key[1] are the read group and qual, respectively
        // and globalDeltaQ is a constant for the read group.  So technically the delta Q reported is simply a lookup
        // into a matrix indexed by rgGroup, qual, and event type.
        // the code below actually creates this cache with a NestedIntegerArray calling into the actual
        // calculateDeltaQReported code.
        final NestedIntegerArray<RecalDatum> byQualTable = recalibrationTables.getQualityScoreTable();
        deltaQReporteds = new NestedIntegerArray<Double>( byQualTable.getDimensions() );
        logger.info("Calculating delta Q reported table...");
        for ( NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int qual = leaf.keys[1];
            final int eventIndex = leaf.keys[2];
            final EventType event = EventType.eventFrom(eventIndex);
            final double globalDeltaQ = getGlobalDeltaQ(rgKey, event);
            final double deltaQReported = calculateDeltaQReported(rgKey, qual, event, globalDeltaQ, (byte)qual);
            deltaQReporteds.put(deltaQReported, rgKey, qual, eventIndex);
        }

        logger.info("done calculating cache");
    }

    /**
     * Recalibrates the base qualities of a read
     *
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     *
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
     * @param read the read to recalibrate
     */
    public void recalibrateRead(final GATKSAMRecord read) {
        if (emitOriginalQuals && read.getAttribute(SAMTag.OQ.name()) == null) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
            } catch (IllegalArgumentException e) {
                throw new UserException.MalformedBAM(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        final ReadCovariates readCovariates = RecalUtils.computeCovariates(read, requestedCovariates);
        final int readLength = read.getReadLength();

        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
            if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
                read.setBaseQualities(null, errorModel);
                continue;
            }

            final byte[] quals = read.getBaseQualities(errorModel);

            // get the keyset for this base using the error model
            final int[][] fullReadKeySet = readCovariates.getKeySet(errorModel);

            // the rg key is constant over the whole read, the global deltaQ is too
            final int rgKey = fullReadKeySet[0][0];

            final double globalDeltaQ = getGlobalDeltaQ(rgKey, errorModel);

            for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read
                final byte origQual = quals[offset];

                // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                if ( origQual >= preserveQLessThan ) {
                    // get the keyset for this base using the error model
                    final int[] keySet = fullReadKeySet[offset];
                    final double deltaQReported = getDeltaQReported(keySet[0], keySet[1], errorModel, globalDeltaQ);
                    final double deltaQCovariates = calculateDeltaQCovariates(recalibrationTables, keySet, errorModel, globalDeltaQ, deltaQReported, origQual);

                    // calculate the recalibrated qual using the BQSR formula
                    double recalibratedQualDouble = origQual + globalDeltaQ + deltaQReported + deltaQCovariates;

                    // recalibrated quality is bound between 1 and MAX_QUAL
                    final byte recalibratedQual = QualityUtils.boundQual(MathUtils.fastRound(recalibratedQualDouble), QualityUtils.MAX_RECALIBRATED_Q_SCORE);

                    // return the quantized version of the recalibrated quality
                    final byte recalibratedQualityScore = quantizationInfo.getQuantizedQuals().get(recalibratedQual);

                    quals[offset] = recalibratedQualityScore;
                }
            }

            // finally update the base qualities in the read
            read.setBaseQualities(quals, errorModel);
        }
    }

    private double getGlobalDeltaQ(final int rgKey, final EventType errorModel) {
        final Double cached = globalDeltaQs.get(rgKey, errorModel.ordinal());

        if ( TEST_CACHING ) {
            final double calcd = calculateGlobalDeltaQ(rgKey, errorModel);
            if ( calcd != cached )
                throw new IllegalStateException("calculated " + calcd + " and cached " + cached + " global delta q not equal at " + rgKey + " / " + errorModel);
        }

        return cachedWithDefault(cached);
    }

    private double getDeltaQReported(final int rgKey, final int qualKey, final EventType errorModel, final double globalDeltaQ) {
        final Double cached = deltaQReporteds.get(rgKey, qualKey, errorModel.ordinal());

        if ( TEST_CACHING ) {
            final double calcd = calculateDeltaQReported(rgKey, qualKey, errorModel, globalDeltaQ, (byte)qualKey);
            if ( calcd != cached )
                throw new IllegalStateException("calculated " + calcd + " and cached " + cached + " global delta q not equal at " + rgKey + " / " + qualKey + " / " + errorModel);
        }

        return cachedWithDefault(cached);
    }

    /**
     * @param d a Double (that may be null) that is the result of a delta Q calculation
     * @return a double == d if d != null, or 0.0 if it is
     */
    private double cachedWithDefault(final Double d) {
        return d == null ? 0.0 : d;
    }

    /**
     * Note that this calculation is a constant for each rgKey and errorModel.  We need only
     * compute this value once for all data.
     *
     * @param rgKey        read group key
     * @param errorModel   the event type
     * @return global delta Q
     */
    private double calculateGlobalDeltaQ(final int rgKey, final EventType errorModel) {
        double result = 0.0;

        final RecalDatum empiricalQualRG = recalibrationTables.getReadGroupTable().get(rgKey, errorModel.ordinal());

        if (empiricalQualRG != null) {
            final double globalDeltaQEmpirical = empiricalQualRG.getEmpiricalQuality();
            final double aggregrateQReported = empiricalQualRG.getEstimatedQReported();
            result = globalDeltaQEmpirical - aggregrateQReported;
        }

        return result;
    }

    private double calculateDeltaQReported(final int rgKey, final int qualKey, final EventType errorModel, final double globalDeltaQ, final byte qualFromRead) {
        double result = 0.0;

        final RecalDatum empiricalQualQS = recalibrationTables.getQualityScoreTable().get(rgKey, qualKey, errorModel.ordinal());
        if (empiricalQualQS != null) {
            final double deltaQReportedEmpirical = empiricalQualQS.getEmpiricalQuality();
            result = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        return result;
    }

    private double calculateDeltaQCovariates(final RecalibrationTables recalibrationTables, final int[] key, final EventType errorModel, final double globalDeltaQ, final double deltaQReported, final byte qualFromRead) {
        double result = 0.0;

        // for all optional covariates
        for (int i = 2; i < requestedCovariates.length; i++) {
            if (key[i] < 0)
                continue;

            result += calculateDeltaQCovariate(recalibrationTables.getTable(i),
                    key[0], key[1], key[i], errorModel,
                    globalDeltaQ, deltaQReported, qualFromRead);
        }

        return result;
    }

    private double calculateDeltaQCovariate(final NestedIntegerArray<RecalDatum> table,
                                            final int rgKey,
                                            final int qualKey,
                                            final int tableKey,
                                            final EventType errorModel,
                                            final double globalDeltaQ,
                                            final double deltaQReported,
                                            final byte qualFromRead) {
        final RecalDatum empiricalQualCO = table.get(rgKey, qualKey, tableKey, errorModel.ordinal());
        if (empiricalQualCO != null) {
            final double deltaQCovariateEmpirical = empiricalQualCO.getEmpiricalQuality();
            return deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported);
        } else {
            return 0.0;
        }
    }
}

/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import com.google.java.contract.Ensures;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.recalibration.EventType;
import org.broadinstitute.gatk.engine.recalibration.covariates.Covariate;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Collections;
import java.util.Iterator;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
    private final double globalQScorePrior;
    private final boolean emitOriginalQuals;

    private byte[] staticQuantizedMapping = null;

    /**
     * Constructor using a GATK Report file
     *
     * @param RECAL_FILE         a GATK Report file containing the recalibration information
     * @param quantizationLevels number of bins to quantize the quality scores
     * @param disableIndelQuals  if true, do not emit base indel qualities
     * @param preserveQLessThan  preserve quality scores less than this value
     * @param staticQuantizedQuals static quantized bins for quality scores
     */
    public BaseRecalibration(final File RECAL_FILE, final int quantizationLevels, final boolean disableIndelQuals, final int preserveQLessThan, final boolean emitOriginalQuals, final double globalQScorePrior, final List<Integer> staticQuantizedQuals, final boolean roundDown) {
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
        this.globalQScorePrior = globalQScorePrior;
        this.emitOriginalQuals = emitOriginalQuals;

        // staticQuantizedQuals is entirely separate from the dynamic binning that quantizationLevels, and
        // staticQuantizedQuals does not make use of quantizationInfo
        if(staticQuantizedQuals != null) {
            if(staticQuantizedQuals.isEmpty()) {
                throw new IllegalStateException("List of static quantized quals is empty.");
            }
            staticQuantizedMapping = constructStaticQuantizedMapping(staticQuantizedQuals, roundDown);
        }
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
            final RecalDatum empiricalQualRG = recalibrationTables.getReadGroupTable().get(rgKey, errorModel.ordinal());

            if( empiricalQualRG != null ) {
                final double epsilon = ( globalQScorePrior > 0.0 && errorModel.equals(EventType.BASE_SUBSTITUTION) ? globalQScorePrior : empiricalQualRG.getEstimatedQReported() );

                for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read
                    final byte origQual = quals[offset];

                    // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                    if ( origQual >= preserveQLessThan ) {
                        // get the keyset for this base using the error model
                        final int[] keySet = fullReadKeySet[offset];
                        final RecalDatum empiricalQualQS = recalibrationTables.getQualityScoreTable().get(keySet[0], keySet[1], errorModel.ordinal());
                        final List<RecalDatum> empiricalQualCovs = new ArrayList<RecalDatum>();
                        for (int i = 2; i < requestedCovariates.length; i++) {
                            if (keySet[i] < 0) {
                                continue;
                            }
                            empiricalQualCovs.add(recalibrationTables.getTable(i).get(keySet[0], keySet[1], keySet[i], errorModel.ordinal()));
                        }

                        double recalibratedQualDouble = hierarchicalBayesianQualityEstimate( epsilon, empiricalQualRG, empiricalQualQS, empiricalQualCovs );

                        // recalibrated quality is bound between 1 and MAX_QUAL
                        final byte recalibratedQual = QualityUtils.boundQual(MathUtils.fastRound(recalibratedQualDouble), RecalDatum.MAX_RECALIBRATED_Q_SCORE);

                        // return the quantized version of the recalibrated quality
                        final byte recalibratedQualityScore = quantizationInfo.getQuantizedQuals().get(recalibratedQual);

                        // Bin to static quals
                        if(staticQuantizedMapping != null) {
                            quals[offset] = staticQuantizedMapping[recalibratedQualityScore];
                        }
                        else {
                            quals[offset] = recalibratedQualityScore;
                        }
                    }
                }
            }

            // finally update the base qualities in the read
            read.setBaseQualities(quals, errorModel);
        }
    }

    /**
     * Constructs an array that maps particular quantized values to a rounded value in staticQuantizedQuals
     *
     * Rounding is done in probability space.  When roundDown is true, we simply round down to the nearest
     * available qual in staticQuantizedQuals
     *
     * @param staticQuantizedQuals the list of qualities to round to
     * @param roundDown round down if true, round to nearest (in probability space) otherwise
     * @return  Array where index representing the quality score to be mapped and the value is the rounded quality score
     */
    protected static byte[] constructStaticQuantizedMapping(List<Integer> staticQuantizedQuals, boolean roundDown) {
        // Create array mapping that maps quals to their rounded value.
        byte[] mapping = new byte[QualityUtils.MAX_QUAL];

        Collections.sort(staticQuantizedQuals);
        Iterator<Integer> quantizationIterator = staticQuantizedQuals.iterator();

        // Fill mapping with one-to-one mappings for values between 0 and MIN_USABLE_Q_SCORE
        // This ensures that quals used as special codes will be preserved
        for(int i = 0 ; i < QualityUtils.MIN_USABLE_Q_SCORE ; i++) {
            mapping[i] = (byte) i;
        }

        // If only one staticQuantizedQual is given, fill mappings larger than QualityUtils.MAX_QUAL with that value
        if(staticQuantizedQuals.size() == 1) {
            int onlyQual = quantizationIterator.next();
            for(int i = QualityUtils.MIN_USABLE_Q_SCORE ; i < QualityUtils.MAX_QUAL ; i++) {
                mapping[i] = (byte) onlyQual;
            }
            return mapping;
        }

        int firstQual = QualityUtils.MIN_USABLE_Q_SCORE;
        int previousQual = firstQual;
        double previousProb = QualityUtils.qualToProb(previousQual);
        while(quantizationIterator.hasNext()) {
            final int nextQual = quantizationIterator.next();
            final double nextProb = QualityUtils.qualToProb(nextQual);

            for (int i = previousQual ; i < nextQual ; i++) {
                if (roundDown) {
                    mapping[i] = (byte) previousQual;
                } else {
                    final double iProb = QualityUtils.qualToProb(i);
                    if ((iProb - previousProb) > (nextProb - iProb)) {
                        mapping[i] = (byte) nextQual;
                    } else {
                        mapping[i] = (byte) previousQual;
                    }
                }
            }
            previousQual = nextQual;
            previousProb = nextProb;
        }
        // Round all quals larger than the largest static qual down to the largest static qual
        for(int j = previousQual ; j < QualityUtils.MAX_QUAL ; j++) {
            mapping[j] = (byte) previousQual;
        }
        return mapping;
    }

    @Ensures("result > 0.0")
    protected static double hierarchicalBayesianQualityEstimate( final double epsilon, final RecalDatum empiricalQualRG, final RecalDatum empiricalQualQS, final List<RecalDatum> empiricalQualCovs ) {
        final double globalDeltaQ = ( empiricalQualRG == null ? 0.0 : empiricalQualRG.getEmpiricalQuality(epsilon) - epsilon );
        final double deltaQReported = ( empiricalQualQS == null ? 0.0 : empiricalQualQS.getEmpiricalQuality(globalDeltaQ + epsilon) - (globalDeltaQ + epsilon) );
        double deltaQCovariates = 0.0;
        for( final RecalDatum empiricalQualCov : empiricalQualCovs ) {
            deltaQCovariates += ( empiricalQualCov == null ? 0.0 : empiricalQualCov.getEmpiricalQuality(deltaQReported + globalDeltaQ + epsilon) - (deltaQReported + globalDeltaQ + epsilon) );
        }

        return epsilon + globalDeltaQ + deltaQReported + deltaQCovariates;
    }
}

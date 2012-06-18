package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.recalibration.QualQuantizer;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Class that encapsulates the information necessary for quality score quantization for BQSR
 *
 * @author carneiro
 * @since 3/26/12
 */
public class QuantizationInfo {
    private List<Byte> quantizedQuals;
    private List<Long> empiricalQualCounts;
    private int quantizationLevels;

    private QuantizationInfo(List<Byte> quantizedQuals, List<Long> empiricalQualCounts, int quantizationLevels) {
        this.quantizedQuals = quantizedQuals;
        this.empiricalQualCounts = empiricalQualCounts;
        this.quantizationLevels = quantizationLevels;
    }

    public QuantizationInfo(List<Byte> quantizedQuals, List<Long> empiricalQualCounts) {
        this(quantizedQuals, empiricalQualCounts, calculateQuantizationLevels(quantizedQuals));
    }
    
    public QuantizationInfo(Map<BQSRKeyManager, Map<Long, RecalDatum>> keysAndTablesMap, int quantizationLevels) {
        final Long [] qualHistogram = new Long[QualityUtils.MAX_QUAL_SCORE+1];                                          // create a histogram with the empirical quality distribution
        for (int i = 0; i < qualHistogram.length; i++)
            qualHistogram[i] = 0L;

        Map<Long, RecalDatum> qualTable = null;                                                                         // look for the quality score table
        for (Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> entry : keysAndTablesMap.entrySet()) {
            BQSRKeyManager keyManager = entry.getKey();
            if (keyManager.getNumRequiredCovariates() == 2)                                                             // it should be the only one with 2 required covariates
                qualTable = entry.getValue();
        }

        if (qualTable == null)
            throw new ReviewedStingException("Could not find QualityScore table.");

        for (RecalDatum datum : qualTable.values()) {
            int empiricalQual = (int) Math.round(datum.getEmpiricalQuality());                                          // convert the empirical quality to an integer ( it is already capped by MAX_QUAL )
            long nObservations = datum.numObservations;
            qualHistogram[empiricalQual] += nObservations;                                                              // add the number of observations for every key
        }
        empiricalQualCounts = Arrays.asList(qualHistogram);                                                             // histogram with the number of observations of the empirical qualities
        quantizeQualityScores(quantizationLevels);

        this.quantizationLevels = quantizationLevels;
    }


    public void quantizeQualityScores(int nLevels) {
        QualQuantizer quantizer = new QualQuantizer(empiricalQualCounts, nLevels, QualityUtils.MIN_USABLE_Q_SCORE);     // quantize the qualities to the desired number of levels
        quantizedQuals = quantizer.getOriginalToQuantizedMap();                                                         // map with the original to quantized qual map (using the standard number of levels in the RAC)
    }

    public void noQuantization() {
        this.quantizationLevels = QualityUtils.MAX_QUAL_SCORE;
        for (int i = 0; i < this.quantizationLevels; i++)
            quantizedQuals.set(i, (byte) i);
    }

    public List<Byte> getQuantizedQuals() {
        return quantizedQuals;
    }

    public int getQuantizationLevels() {
        return quantizationLevels;
    }

    public GATKReportTable generateReportTable() {
        GATKReportTable quantizedTable = new GATKReportTable(RecalDataManager.QUANTIZED_REPORT_TABLE_TITLE, "Quality quantization map", 3);
        quantizedTable.addColumn(RecalDataManager.QUALITY_SCORE_COLUMN_NAME);
        quantizedTable.addColumn(RecalDataManager.QUANTIZED_COUNT_COLUMN_NAME);
        quantizedTable.addColumn(RecalDataManager.QUANTIZED_VALUE_COLUMN_NAME);

        for (int qual = 0; qual <= QualityUtils.MAX_QUAL_SCORE; qual++) {
            quantizedTable.set(qual, RecalDataManager.QUALITY_SCORE_COLUMN_NAME, qual);
            quantizedTable.set(qual, RecalDataManager.QUANTIZED_COUNT_COLUMN_NAME, empiricalQualCounts.get(qual));
            quantizedTable.set(qual, RecalDataManager.QUANTIZED_VALUE_COLUMN_NAME, quantizedQuals.get(qual));
        }
        return quantizedTable;
    }

    private static int calculateQuantizationLevels(List<Byte> quantizedQuals) {
        byte lastByte = -1;
        int quantizationLevels = 0;
        for (byte q : quantizedQuals) {
            if (q != lastByte) {
                quantizationLevels++;
                lastByte = q;
            }
        }
        return quantizationLevels;
    }
}

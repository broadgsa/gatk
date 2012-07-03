package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for on-the-fly recalibration.
 *
 * @author Mauricio Carneiro
 * @since 3/16/12
 */
public class BaseRecalibrationUnitTest {

    private org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager dataManager;

    private ReadGroupCovariate rgCovariate;
    private QualityScoreCovariate qsCovariate;
    private ContextCovariate cxCovariate;
    private CycleCovariate cyCovariate;

    private GATKSAMRecord read = ReadUtils.createRandomRead(10000);
    private BaseRecalibration baseRecalibration;
    private ReadCovariates readCovariates;


    @BeforeClass
    public void init() {
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("rg");
        rg.setPlatform("illumina");
        read.setReadGroup(rg);

        byte[] quals = new byte[read.getReadLength()];
        for (int i = 0; i < read.getReadLength(); i++)
            quals[i] = 20;
        read.setBaseQualities(quals);

        RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        List<Covariate> requiredCovariates = new ArrayList<Covariate>();
        List<Covariate> optionalCovariates = new ArrayList<Covariate>();

        dataManager = new org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager(true, 4);

        rgCovariate = new ReadGroupCovariate();
        rgCovariate.initialize(RAC);
        requiredCovariates.add(rgCovariate);

        qsCovariate = new QualityScoreCovariate();
        qsCovariate.initialize(RAC);
        requiredCovariates.add(qsCovariate);

        cxCovariate = new ContextCovariate();
        cxCovariate.initialize(RAC);
        optionalCovariates.add(cxCovariate);
        cyCovariate = new CycleCovariate();
        cyCovariate.initialize(RAC);
        optionalCovariates.add(cyCovariate);

        final Covariate[] requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate cov : requiredCovariates)
            requestedCovariates[covariateIndex++] = cov;
        for (final Covariate cov : optionalCovariates)
            requestedCovariates[covariateIndex++] = cov;

        readCovariates = RecalDataManager.computeCovariates(read, requestedCovariates);

        RecalibrationTables recalibrationTables = new RecalibrationTables(requestedCovariates);
        final NestedIntegerArray<RecalDatum> rgTable = recalibrationTables.getTable(RecalibrationTables.TableType.READ_GROUP_TABLE);
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getTable(RecalibrationTables.TableType.QUALITY_SCORE_TABLE);

        for (int i=0; i<read.getReadLength(); i++) {
            final int[] bitKeys = readCovariates.getMismatchesKeySet(i);
            final Object[] objKey = buildObjectKey(bitKeys);

            Random random = new Random();
            int nObservations = random.nextInt(10000);
            int nErrors = random.nextInt(10);
            double estimatedQReported = 30;
            double empiricalQuality = calcEmpiricalQual(nObservations, nErrors);

            org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum oldDatum = new org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);
            dataManager.addToAllTables(objKey, oldDatum, QualityUtils.MIN_USABLE_Q_SCORE);

            RecalDatum newDatum = new RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);

            rgTable.put(newDatum, bitKeys[0], EventType.BASE_SUBSTITUTION.index);
            qualTable.put(newDatum, bitKeys[0], bitKeys[1], EventType.BASE_SUBSTITUTION.index);
            for (int j = 0; j < optionalCovariates.size(); j++) {
                final NestedIntegerArray<RecalDatum> covTable = recalibrationTables.getTable(RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.index + j);
                covTable.put(newDatum, bitKeys[0], bitKeys[1], j, bitKeys[RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.index + j], EventType.BASE_SUBSTITUTION.index);
            }
        }

    dataManager.generateEmpiricalQualities(1, QualityUtils.MAX_RECALIBRATED_Q_SCORE);

        List<Byte> quantizedQuals = new ArrayList<Byte>();
        List<Long> qualCounts = new ArrayList<Long>();
        for (byte i = 0; i <= QualityUtils.MAX_QUAL_SCORE; i++) {
            quantizedQuals.add(i);
            qualCounts.add(1L);
        }
        QuantizationInfo quantizationInfo = new QuantizationInfo(quantizedQuals, qualCounts);
        quantizationInfo.noQuantization();
        baseRecalibration = new BaseRecalibration(quantizationInfo, recalibrationTables, requestedCovariates);

    }


    @Test(enabled=false)
    public void testGoldStandardComparison() {
        for (int i = 0; i < read.getReadLength(); i++) {
            int [] bitKey = readCovariates.getKeySet(i, EventType.BASE_SUBSTITUTION);
            Object [] objKey = buildObjectKey(bitKey);
            byte v2 = baseRecalibration.performSequentialQualityCalculation(bitKey, EventType.BASE_SUBSTITUTION);
            byte v1 = goldStandardSequentialCalculation(objKey);
            Assert.assertEquals(v2, v1);
        }
    }

    private Object[] buildObjectKey(final int[] bitKey) {
        Object[] key = new Object[bitKey.length];
        key[0] = rgCovariate.formatKey(bitKey[0]);
        key[1] = qsCovariate.formatKey(bitKey[1]);
        key[2] = cxCovariate.formatKey(bitKey[2]);
        key[3] = cyCovariate.formatKey(bitKey[3]);
        return key;
    }

    private static void printNestedHashMap(Map table, String output) {
        for (Object key : table.keySet()) {
            String ret;
            if (output.isEmpty())
                ret = "" + key;
            else
                ret = output + "," + key;

            Object next = table.get(key);
            if (next instanceof org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum)
                System.out.println(ret + " => " + next);
            else
                printNestedHashMap((Map) next, "" + ret);
        }
    }

    private void updateCovariateWithKeySet(final Map<Long, RecalDatum> recalTable, final Long hashKey, final RecalDatum datum) {
        RecalDatum previousDatum = recalTable.get(hashKey);                                                             // using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        if (previousDatum == null)                                                                                      // key doesn't exist yet in the map so make a new bucket and add it
            recalTable.put(hashKey, datum.copy());
        else
            previousDatum.combine(datum);                                                                               // add one to the number of observations and potentially one to the number of mismatches
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
     * @param key The list of Comparables that were calculated from the covariates
     * @return A recalibrated quality score as a byte
     */
    private byte goldStandardSequentialCalculation(final Object... key) {

        final byte qualFromRead = (byte) Integer.parseInt(key[1].toString());
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] qualityScoreCollapsedKey = new Object[2];
        final Object[] covariateCollapsedKey = new Object[3];

        // The global quality shift (over the read group only)
        readGroupCollapsedKey[0] = key[0];
        final org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum globalRecalDatum = ((org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum) dataManager.getCollapsedTable(0).get(readGroupCollapsedKey));
        double globalDeltaQ = 0.0;
        if (globalRecalDatum != null) {
            final double globalDeltaQEmpirical = globalRecalDatum.getEmpiricalQuality();
            final double aggregrateQReported = globalRecalDatum.getEstimatedQReported();
            globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
        }

        // The shift in quality between reported and empirical
        qualityScoreCollapsedKey[0] = key[0];
        qualityScoreCollapsedKey[1] = key[1];
        final org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum qReportedRecalDatum = ((org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum) dataManager.getCollapsedTable(1).get(qualityScoreCollapsedKey));
        double deltaQReported = 0.0;
        if (qReportedRecalDatum != null) {
            final double deltaQReportedEmpirical = qReportedRecalDatum.getEmpiricalQuality();
            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        // The shift in quality due to each covariate by itself in turn
        double deltaQCovariates = 0.0;
        double deltaQCovariateEmpirical;
        covariateCollapsedKey[0] = key[0];
        covariateCollapsedKey[1] = key[1];
        for (int iii = 2; iii < key.length; iii++) {
            covariateCollapsedKey[2] = key[iii]; // The given covariate
            final org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum covariateRecalDatum = ((org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum) dataManager.getCollapsedTable(iii).get(covariateCollapsedKey));
            if (covariateRecalDatum != null) {
                deltaQCovariateEmpirical = covariateRecalDatum.getEmpiricalQuality();
                deltaQCovariates += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
            }
        }

        final double newQuality = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
        return QualityUtils.boundQual((int) Math.round(newQuality), QualityUtils.MAX_RECALIBRATED_Q_SCORE);

        // Verbose printouts used to validate with old recalibrator
        //if(key.contains(null)) {
        //    System.out.println( key  + String.format(" => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte));
        //}
        //else {
        //    System.out.println( String.format("%s %s %s %s => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 key.get(0).toString(), key.get(3).toString(), key.get(2).toString(), key.get(1).toString(), qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte) );
        //}

        //return newQualityByte;
    }

    public static double calcEmpiricalQual(final int observations, final int errors) {
        final int smoothing = 1;
        final double doubleMismatches = (double) (errors + smoothing);
        final double doubleObservations = (double) ( observations + smoothing );
        double empiricalQual = -10 * Math.log10(doubleMismatches / doubleObservations);
        return Math.min(QualityUtils.MAX_RECALIBRATED_Q_SCORE, empiricalQual);
    }
}

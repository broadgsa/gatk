package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
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
    private LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>> keysAndTablesMap;

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
        ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>();

        dataManager = new org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager(true, 4);
        keysAndTablesMap = new LinkedHashMap<BQSRKeyManager, Map<Long, RecalDatum>>();

        rgCovariate = new ReadGroupCovariate();
        rgCovariate.initialize(RAC);
        requiredCovariates.add(rgCovariate);
        BQSRKeyManager rgKeyManager = new BQSRKeyManager(requiredCovariates, optionalCovariates);
        keysAndTablesMap.put(rgKeyManager, new HashMap<Long, RecalDatum>());

        qsCovariate = new QualityScoreCovariate();
        qsCovariate.initialize(RAC);
        requiredCovariates.add(qsCovariate);
        BQSRKeyManager qsKeyManager = new BQSRKeyManager(requiredCovariates, optionalCovariates);
        keysAndTablesMap.put(qsKeyManager, new HashMap<Long, RecalDatum>());

        cxCovariate = new ContextCovariate();
        cxCovariate.initialize(RAC);
        optionalCovariates.add(cxCovariate);
        cyCovariate = new CycleCovariate();
        cyCovariate.initialize(RAC);
        optionalCovariates.add(cyCovariate);
        BQSRKeyManager cvKeyManager = new BQSRKeyManager(requiredCovariates, optionalCovariates);
        keysAndTablesMap.put(cvKeyManager, new HashMap<Long, RecalDatum>());

        for (Covariate cov : requiredCovariates)
            requestedCovariates.add(cov);
        for (Covariate cov : optionalCovariates)
            requestedCovariates.add(cov);

        readCovariates = RecalDataManager.computeCovariates(read, requestedCovariates);

        for (int i=0; i<read.getReadLength(); i++) {
            Long[] bitKeys = readCovariates.getMismatchesKeySet(i);

            Object[] objKey = buildObjectKey(bitKeys);

            Random random = new Random();
            int nObservations = random.nextInt(10000);
            int nErrors = random.nextInt(10);
            double estimatedQReported = 30;
            double empiricalQuality = calcEmpiricalQual(nObservations, nErrors);

            org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum oldDatum = new org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);
            dataManager.addToAllTables(objKey, oldDatum, QualityUtils.MIN_USABLE_Q_SCORE);

            RecalDatum newDatum = new RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);
            for (Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> mapEntry : keysAndTablesMap.entrySet()) {
                Long[] keys = mapEntry.getKey().longsFromAllKeys(bitKeys, EventType.BASE_SUBSTITUTION);
                for (Long key : keys)
                    updateCovariateWithKeySet(mapEntry.getValue(), key, newDatum);
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
        baseRecalibration = new BaseRecalibration(quantizationInfo, keysAndTablesMap, requestedCovariates.toArray());

    }


    @Test(enabled=false)
    public void testGoldStandardComparison() {
        debugTables();
        for (int i = 0; i < read.getReadLength(); i++) {
            Long [] bitKey = readCovariates.getKeySet(i, EventType.BASE_SUBSTITUTION);
            Object [] objKey = buildObjectKey(bitKey);
            byte v2 = baseRecalibration.performSequentialQualityCalculation(bitKey, EventType.BASE_SUBSTITUTION);
            byte v1 = goldStandardSequentialCalculation(objKey);
            Assert.assertEquals(v2, v1);
        }
    }

    private Object[] buildObjectKey(Long[] bitKey) {
        Object[] key = new Object[bitKey.length];
        key[0] = rgCovariate.formatKey(bitKey[0]);
        key[1] = qsCovariate.formatKey(bitKey[1]);
        key[2] = cxCovariate.formatKey(bitKey[2]);
        key[3] = cyCovariate.formatKey(bitKey[3]);
        return key;
    }

    private void debugTables() {
        System.out.println("\nV1 Table\n");
        System.out.println("ReadGroup Table:");
        NestedHashMap nestedTable = dataManager.getCollapsedTable(0);
        printNestedHashMap(nestedTable.data, "");
        System.out.println("\nQualityScore Table:");
        nestedTable = dataManager.getCollapsedTable(1);
        printNestedHashMap(nestedTable.data, "");
        System.out.println("\nCovariates Table:");
        nestedTable = dataManager.getCollapsedTable(2);
        printNestedHashMap(nestedTable.data, "");
        nestedTable = dataManager.getCollapsedTable(3);
        printNestedHashMap(nestedTable.data, "");


        int i = 0;
        System.out.println("\nV2 Table\n");
        for (Map.Entry<BQSRKeyManager, Map<Long, RecalDatum>> mapEntry : keysAndTablesMap.entrySet()) {
            BQSRKeyManager keyManager = mapEntry.getKey();
            Map<Long, RecalDatum> table = mapEntry.getValue();
            switch(i++) {
                case 0 :
                    System.out.println("ReadGroup Table:");
                    break;
                case 1 :
                    System.out.println("QualityScore Table:");
                    break;
                case 2 :
                    System.out.println("Covariates Table:");
                    break;
            }
            for (Map.Entry<Long, RecalDatum> entry : table.entrySet()) {
                Long key = entry.getKey();
                RecalDatum datum = entry.getValue();
                List<Object> keySet = keyManager.keySetFrom(key);
                System.out.println(String.format("%s => %s", Utils.join(",", keySet), datum) + "," + datum.getEstimatedQReported());
            }
            System.out.println();
        }


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

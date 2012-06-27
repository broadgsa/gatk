package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author carneiro
 * @since 4/21/12
 */
public class RecalibrationReportUnitTest {
    @Test(enabled = false)
    public void testOutput() {
        final int length = 100;

        List<Byte> quals = new ArrayList<Byte>(QualityUtils.MAX_QUAL_SCORE + 1);
        List<Long> counts = new ArrayList<Long>(QualityUtils.MAX_QUAL_SCORE + 1);

        for (int i = 0;  i<= QualityUtils.MAX_QUAL_SCORE; i++) {
            quals.add((byte) i);
            counts.add(1L);
        }

        final QuantizationInfo quantizationInfo = new QuantizationInfo(quals, counts);
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        quantizationInfo.noQuantization();
        final List<Covariate> requiredCovariates = new LinkedList<Covariate>();
        final List<Covariate> optionalCovariates = new LinkedList<Covariate>();

        final ReadGroupCovariate rgCovariate = new ReadGroupCovariate();
        rgCovariate.initialize(RAC);
        requiredCovariates.add(rgCovariate);

        final QualityScoreCovariate qsCovariate = new QualityScoreCovariate();
        qsCovariate.initialize(RAC);
        requiredCovariates.add(qsCovariate);

        final ContextCovariate cxCovariate = new ContextCovariate();
        cxCovariate.initialize(RAC);
        optionalCovariates.add(cxCovariate);
        final CycleCovariate cyCovariate = new CycleCovariate();
        cyCovariate.initialize(RAC);
        optionalCovariates.add(cyCovariate);

        final Covariate[] requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate cov : requiredCovariates)
            requestedCovariates[covariateIndex++] = cov;
        for (final Covariate cov : optionalCovariates)
            requestedCovariates[covariateIndex++] = cov;

        final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("id");
        rg.setPlatform("illumina");
        final GATKSAMRecord read = ReadUtils.createRandomRead(length, false);
        read.setReadGroup(rg);
        final byte [] readQuals = new byte[length];
        for (int i = 0; i < length; i++)
            readQuals[i] = 20;
        read.setBaseQualities(readQuals);

        final int expectedKeys = expectedNumberOfKeys(4, length, RAC.INDELS_CONTEXT_SIZE, RAC.MISMATCHES_CONTEXT_SIZE);
        int nKeys = 0;                                                                                                  // keep track of how many keys were produced
        final ReadCovariates rc = RecalDataManager.computeCovariates(read, requestedCovariates);

        final NestedHashMap rgTable = new NestedHashMap();
        final NestedHashMap qualTable = new NestedHashMap();
        final NestedHashMap covTable = new NestedHashMap();

        for (int offset = 0; offset < length; offset++) {

            for (EventType errorMode : EventType.values()) {

                final int[] covariates = rc.getKeySet(offset, errorMode);
                final int randomMax = errorMode == EventType.BASE_SUBSTITUTION ? 10000 : 100000;

                rgTable.put(RecalDatum.createRandomRecalDatum(randomMax, 10), covariates[0], errorMode.index);
                qualTable.put(RecalDatum.createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], errorMode.index);
                nKeys += 2;
                for (int j = 0; j < optionalCovariates.size(); j++) {
                    covTable.put(RecalDatum.createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], j, covariates[2 + j], errorMode.index);
                    nKeys++;
                }
            }
        }
        Assert.assertEquals(nKeys, expectedKeys);

        final RecalibrationTables recalibrationTables = new RecalibrationTables(rgTable, qualTable, covTable);

        final RecalibrationReport report = new RecalibrationReport(quantizationInfo, recalibrationTables, RAC.generateReportTable(), RAC);

        File output = new File("RecalibrationReportUnitTestOutuput.grp");
        PrintStream out;
        try {
            out = new PrintStream(output);
        } catch (FileNotFoundException e) {
            throw new ReviewedStingException("couldn't create the file " + output, e);
        }
        report.output(out);

        RecalibrationReport loadedReport = new RecalibrationReport(output);

        Assert.assertTrue(report.equals(loadedReport));
        if (!output.delete())
            throw new ReviewedStingException("File could not be deleted " + output);
    }

    private static int expectedNumberOfKeys (int nCovariates, int readLength, int indelContextSize, int mismatchesContextSize) {
        int nommcs = readLength >= mismatchesContextSize ? mismatchesContextSize-1 : readLength;
        int noincs = readLength >= indelContextSize ? 2*(indelContextSize-1) : 2*readLength;
        return (nCovariates * readLength  * 3) -  nommcs - noincs;
    }

}

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 8/20/12
 * Time: 9:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class VariantsToBinaryPedIntegrationTest extends WalkerTest {

    public static final String VTBP_DATA_DIR =  "/humgen/gsa-hpprojects/GATK/data/Validation_Data/VariantsToBinaryPed/";

    public static String baseTestString(String inputVCF, String inputMetaData, int gq) {
        return "-T VariantsToBinaryPed -R " + b37KGReference +
                " -V " + VTBP_DATA_DIR+inputVCF + " -m "+VTBP_DATA_DIR+inputMetaData + String.format(" -mgq %d",gq) +
                " -bim %s -fam %s -bed %s";

    }

    @Test
    public void testNA12878Alone() {
        String testName = "testNA12878Alone";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12878.subset.vcf", "CEUTrio.NA12878.fam",10),
                3,
                Arrays.asList("411ef932095728bfa5e509c2c0e4cfa8","8e8bc0b5e69f22c54c0960f13c25d26c","02f1c462ebc8576e399d0e94f729fd95")
        );

        executeTest(testName, spec);
    }

    @Test
    public void testNA12878AloneMetaData() {
        String testName = "testNA12878AloneMetaData";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12878.subset.vcf", "CEUTrio.NA12878.metadata.txt",10),
                3,
                Arrays.asList("411ef932095728bfa5e509c2c0e4cfa8","7251ca4e8a515b698e7e7d25cff91978","02f1c462ebc8576e399d0e94f729fd95")
        );

        executeTest(testName, spec);
    }

    @Test
    public void testCEUTrio() {
        String testName = "testCEUTrio";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("CEUTrio.subset.vcf", "CEUTrio.fam",10),
                3,
                Arrays.asList("59b93fbb4bb31309b3adc83ba96dd1a2","900f22c6d49a6ba0774466e99592e51d","7887d2e0bf605dbcd0688c552cdb99d5")
        );

        executeTest(testName, spec);
    }

    @Test
    public void testCEUTrioMetaData() {
        String testName = "testCEUTrioMetaData";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("CEUTrio.subset.vcf", "CEUTrio.metadata.txt",10),
                3,
                Arrays.asList("59b93fbb4bb31309b3adc83ba96dd1a2","2113d2cc0a059e35b1565196b7c5d98f","7887d2e0bf605dbcd0688c552cdb99d5")
        );

        executeTest(testName, spec);
    }

    @Test
    public void testMalformedFam() {
        String testName = "testMalformedFam";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("CEUTrio.subset.vcf", "CEUTrio.malformed.fam",10),
                3,
                UserException.class
        );

        executeTest(testName, spec);
    }
}



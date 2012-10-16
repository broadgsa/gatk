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

    public static String baseTestString(String inputVCF, String inputMetaData, int gq, String mode) {
        return "-T VariantsToBinaryPed -R " + b37KGReference + " -mode "+mode +
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
    public void testNA12878AloneSNPMajor() {
        String testName = "testNA12878AloneSNPMajor";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12878.subset.vcf", "CEUTrio.NA12878.metadata.txt",10,"SNP_MAJOR"),
                3,
                Arrays.asList("411ef932095728bfa5e509c2c0e4cfa8","7251ca4e8a515b698e7e7d25cff91978","ada1acc475d096012b921b3219c3a446")
        );

        executeTest(testName, spec);
    }

    @Test
    public void testNA12878HighGQ() {
        String testName = "testNA12878HighGQ";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12878.subset.vcf", "CEUTrio.NA12878.metadata.txt",80),
                3,
                Arrays.asList("411ef932095728bfa5e509c2c0e4cfa8","7251ca4e8a515b698e7e7d25cff91978","0822adea688e99bb336afe5172d4c959")
        );

        executeTest(testName, spec);
    }

    @Test
    public void testVCFMismatchReference() {
        String testName = "testVCFMismatchReference";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12878.badReference.vcf", "CEUTrio.NA12878.metadata.txt",80),
                3,
                UserException.class
        );

        executeTest(testName, spec);
    }

    @Test
    public void test1000GWithIndels() {
        String testName = "test1000GWithIndels";
         WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("1000G_selected_allVariants.vcf", "1000G_selected_allVariants.md.txt",0),
                3,
                Arrays.asList("3c98112434d9948dc47da72ad14e8d84","3aceda4f9bb5b5457797c1fe5a85b03d","451498ceff06c1649890900fa994f1af")
        );
    }

    @Test
    public void test1000GWithIndelsSNPMajor() {
        String testName = "test1000GWithIndelsSNPMajor";
         WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("1000G_selected_allVariants.vcf", "1000G_selected_allVariants.md.txt",0,"SNP_MAJOR"),
                3,
                Arrays.asList("3c98112434d9948dc47da72ad14e8d84","4a0ba3d0594b06306aa6459e4e28ec9a","451498ceff06c1649890900fa994f1af")
        );
    }

    @Test
    public void test1000G_Symbolic() {
        String testName = "test1000G_Symbolic";
         WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("1000G_selected_SVs.vcf", "1000G_selected_allVariants.md.txt",0),
                3,
                Arrays.asList("5e7ede48e7c5d5972c59dc5558a06e40","451498ceff06c1649890900fa994f1af","4b53a82a0b2d1a22a6eebca50a4f83a8")
        );
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

    @Test
    public void testFailFast() {
        String testName = "testFailFast";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("HapMap.testFailFast.vcf", "HapMap_only_famids.fam",10),
                3,
                UserException.class
        );

        executeTest(testName, spec);
    }

    @Test
    public void testFailFastMeta() {
    String testName = "testFailFastMeta";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("HapMap.testFailFast.vcf", "HapMap_only_famids.metadata.txt",10),
                3,
                UserException.class
        );

        executeTest(testName, spec);

    }

}



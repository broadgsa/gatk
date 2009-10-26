package org.broadinstitute.sting.playground.gatk.walkers;

import org.junit.Test;
import org.broadinstitute.sting.WalkerTest;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 19, 2009
 * Time: 1:57:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class BaseTransitionTableCalculatorJavaIntegrationTest extends WalkerTest{
    // MD5s last computed 10/26 at revision 1897
    public static final String OUTPUT_MD5_STANDARD = "bb7f7c8b71a2a19dbbd47699708816b0";
    public static final String OUTPUT_MD5_3MISMATCHES = "000bd16dfea8df415e2104fe894aec83";
    public static final String OUTPUT_MD5_LOWMAPPINGQUALITY = "db3af69e61e90c59a3ca0ca25605fa67";
    public static final String OUTPUT_MD5_LOWQUALITYSCORE = "f990bd4ba5a951b16603131906b74326";
    public static final String OUTPUT_MD5_LOWCONFIDENTREFTHRESHOLD = "7c6c2ed2110ee3030fcd060346623fcc";
    public static final String OUTPUT_MD5_HIGHCONFIDENTREFTHRESHOLD = "b30efed02e8f6b356e56e5875db40f2c";
    public static final String OUTPUT_MD5_ALLARGUMENTS = "3f1901a40e79300da4cbab1488c4839d";
    public static final String OUTPUT_MD5_USEREADGROUP = "7422ea018b8dea52398b81502cbe5c38";
    public static final String LOCUS = "1:10,000,000-10,200,000";
    public static final String BAM_FILE = "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam";
    public static final String REFERENCE = "/broad/1KG/reference/human_b36_both.fasta";

    @Test
    public void testBaseTransitionTableCalculatorJava() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE;
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_STANDARD));
        executeTest("testBaseTransitionTableCalculatorJava", spec);
    }

    @Test
    public void testBaseTransitionTableCalculatorJavaTreeReduce() {
        String args = "-T BaseTransitionTableCalculatorJava -of %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+" -nt 5";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_STANDARD));
        executeTest("testBaseTransitionTableCalculatorJava: tree reduce", spec);
    }

    @Test
    public void testBaseTransitionTableCalculatorJavaAdditionalMismatches() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+" --maxNumMismatches 3";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_3MISMATCHES));
        executeTest("BaseTransitionTableCalculatorJava: Additional Mismatches",spec);
    }

    @Test
    public void testBaseTransitionCalculatorJavaLowMapingQuality() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+ " --minMappingQuality 5";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_LOWMAPPINGQUALITY));
        executeTest("BaseTransitionTableCalculatorJava: Low Mapping Quality",spec);
    }

    @Test
    public void testBaseTransitionCalculatorJavaLowQualityScore() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+ " --minQualityScore 5";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_LOWQUALITYSCORE));
        executeTest("BaseTransitionTableCalculatorJava: Low Quality Score",spec);
    }

    @Test
    public void testBaseTransitionCalculatorJavaLowConfidentRefThreshold() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+ " --confidentRefThreshold 2";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_LOWCONFIDENTREFTHRESHOLD));
        executeTest("BaseTransitionTableCalculatorJava: Low Ref Threshold",spec);
    }

    @Test
    public void testBaseTransitionCalculatorJavaHighConfidentRefThreshold() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+ " --confidentRefThreshold 8";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_HIGHCONFIDENTREFTHRESHOLD));
        executeTest("BaseTransitionTableCalculatorJava: High Ref Threshold",spec);
    }

    @Test
    public void testBaseTransitionCalculatorJavaAllArguments() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+ " --confidentRefThreshold 2 --minQualityScore 5 --minMappingQuality 5 --maxNumMismatches 3";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_ALLARGUMENTS));
        executeTest("BaseTransitionTableCalculatorJava: Low Ref Threshold, Low Q Score, Low Mapping Quality, Additional Mismatches",spec);
    }

    @Test
    public void testBaseTransitionTableCalculatorJavaUseReadGroup() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+" --useReadGroup";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_USEREADGROUP));
        executeTest("testBaseTransitionTableCalculatorJava: Use read group", spec);
    }
}

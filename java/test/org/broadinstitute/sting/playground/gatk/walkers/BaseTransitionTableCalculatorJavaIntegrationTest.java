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
    // MD5s last computed 10/20 at revision 1877
    public static final String OUTPUT_MD5_STANDARD = "e88f010ea842bcdb6503a4da24b90edc";
    public static final String OUTPUT_MD5_3MISMATCHES = "46f9aadbfe260a286fb6c8cac137dddd";
    public static final String OUTPUT_MD5_LOWMAPPINGQUALITY = "0b7447e0a271ffa5c8ff1719db3585e6";
    public static final String OUTPUT_MD5_LOWQUALITYSCORE = "87d793f751328b0a9299a69b99cd0112";
    public static final String OUTPUT_MD5_LOWCONFIDENTREFTHRESHOLD = "bc70acf997295a6cc13f3e88b254cc24";
    public static final String OUTPUT_MD5_HIGHCONFIDENTREFTHRESHOLD = "17d49c12fc64926a71bc32ac1eec7f68";
    public static final String OUTPUT_MD5_ALLARGUMENTS = "5599a4a577927fb9875adc26140dee7e";
    public static final String OUTPUT_MD5_USEREADGROUP = "f6e8bd2101316deae93512dd09a69567";
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

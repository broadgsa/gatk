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
    public static final String OUTPUT_MD5_STANDARD = "856a3124fb078b4d820f935369ec53cc";
    public static final String OUTPUT_MD5_3MISMATCHES = "847375e3d40415d237c627a56d163760";
    public static final String OUTPUT_MD5_LOWMAPPINGQUALITY = "d7e490168dfef8eb4185b596f55bfaa8";
    public static final String OUTPUT_MD5_LOWQUALITYSCORE = "e3f3c3875936b0df7a8899c2ef5d1dd2";
    public static final String OUTPUT_MD5_LOWCONFIDENTREFTHRESHOLD = "521a9fa7716ed22550c2ba3fe3409070";
    public static final String OUTPUT_MD5_HIGHCONFIDENTREFTHRESHOLD = "8ab6d389fc494881736e9a58126c2f1b";
    public static final String OUTPUT_MD5_ALLARGUMENTS = "f45481946d7a5c70078d432b0baff083";
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
        executeTest("BaseTransitionTableCalculatorJava: Low Ref Threshold",spec);
    }

    @Test
    public void testBaseTransitionCalculatorJavaAllArguments() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE+ " --confidentRefThreshold 2 --minQualityScore 5 --minMappingQuality 5 --maxNumMismatches 3";
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5_ALLARGUMENTS));
        executeTest("BaseTransitionTableCalculatorJava: Low Ref Threshold, Low Q Score, Low Mapping Quality, Additional Mismatches",spec);
    }
}

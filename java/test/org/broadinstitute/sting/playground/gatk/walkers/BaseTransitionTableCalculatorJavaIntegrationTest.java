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

    public static final String OUTPUT_MD5 = "856a3124fb078b4d820f935369ec53cc";
    public static final String LOCUS = "1:10,000,000-10,200,000";
    public static final String BAM_FILE = "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam";
    public static final String REFERENCE = "/broad/1KG/reference/human_b36_both.fasta";

    @Test
    public void testBaseTransitionTableCalculatorJava() {
        String args = "-T BaseTransitionTableCalculatorJava -o %s -I "+BAM_FILE+" -L "+LOCUS+" -R "+REFERENCE;
        WalkerTest.WalkerTestSpec spec =  new WalkerTest.WalkerTestSpec(args,1,Arrays.asList(OUTPUT_MD5));
        executeTest("testBaseTransitionTableCalculatorJava", spec);
    }
}

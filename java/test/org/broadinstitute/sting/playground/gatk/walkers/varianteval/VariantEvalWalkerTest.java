package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author aaron
 *         <p/>
 *         Class VariantEvalWalkerTest
 *         <p/>
 *         test out the variant eval walker under a bunch of different runtime conditions.
 */
public class VariantEvalWalkerTest extends WalkerTest {

    @Test
    public void emptyTest() {

    }

    @Test
    public void testEvalVariantROD() {
        List <String> md5 = new ArrayList<String>();
        md5.add("094c0adb8e4ae4de424f26482fd43152");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /broad/1KG/reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantEval" +
                        " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                        " -L 1:10,000,000-11,000,000" +
                        " --outerr %s" +
                        " --supressDateInformation",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalVariantROD", spec).getFirst();

    }
}


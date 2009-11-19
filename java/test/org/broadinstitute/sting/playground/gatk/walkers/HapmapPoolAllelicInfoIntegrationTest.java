package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 18, 2009
 * Time: 9:29:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class HapmapPoolAllelicInfoIntegrationTest extends WalkerTest {

    @Test
    public void testFHSPool3() {
        String test_args = "-T HapmapPoolAllelicInfo -samples /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_samples.txt "
                + "-B /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_sample_paths.txt "
                + "-B calls,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_raw_calls.geli "
                + "-I /humgen/gsa-scr1/GATK_Data/Validation_Data/FHSP_pool3_test.bam "
                + "-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -of %s "
                + "-ps 40 -L /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_test_intervals.interval_list";
        String md5ForThisTest = "da1222d7514f247eae066134d8e3cde3";
        WalkerTestSpec spec = new WalkerTestSpec(test_args, 1, Arrays.asList(md5ForThisTest));
        executeTest("Pool 3 of FHS Pilot on testbed intervals", spec);
    }

    @Test
    public void testFHSPool3NoIntervals() {
       String test_args = "-T HapmapPoolAllelicInfo -samples /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_samples.txt "
                + "-B /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_sample_paths.txt "
                + "-B calls,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_raw_calls.geli "
                + "-I /humgen/gsa-scr1/GATK_Data/Validation_Data/FHSP_pool3_test.bam "
                + "-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -of %s "
                + "-ps 40";
        String md5ForThisTest = "120f3307d94d613c3559a1051fe3aaef";
        WalkerTestSpec spec = new WalkerTestSpec(test_args, 1, Arrays.asList(md5ForThisTest));
        executeTest("Pool 3 of FHS Pilot without testbed intervals", spec);
    }

    @Test
    public void testSmallPool() {
        String test_args = "-T HapmapPoolAllelicInfo -samples /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878_NA12891_samples.txt "
                +"-B NA12878,GFF,/humgen/gsa-scr1/andrewk/hapmap_1kg/gffs/NA12878.gff "
                +"-B NA12891,GFF,/humgen/gsa-scr1/andrewk/hapmap_1kg/gffs/NA12891.gff "
                +"-I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12891.a2b.recal.annotation_subset.bam "
                +"-L chr1:14000000-18000000 -B calls,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12891.a2b.ssg1b.geli.calls "
                +"-of %s -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -ps 2";
        String md5ForThisTest = "cc81a9d8279e7b6ad50f9dc1f19a7419";
        WalkerTestSpec spec = new WalkerTestSpec(test_args, 1, Arrays.asList(md5ForThisTest));
        executeTest("Test on pool of two individuals -- CEU father+daughter chips against subset of the father's calls; power at lod 3", spec);

    }
}

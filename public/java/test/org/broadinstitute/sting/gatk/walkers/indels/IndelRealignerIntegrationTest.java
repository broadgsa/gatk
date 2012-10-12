package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class IndelRealignerIntegrationTest extends WalkerTest {

    private static final String mainTestBam = validationDataLocation + "indelRealignerTest.pilot1.ceu.fixed.fixmates.bam";
    private static final String mainTestIntervals = validationDataLocation + "indelRealignerTest.pilot1.ceu.intervals";
    private static final String knownIndels = validationDataLocation + "indelRealignerTest.pilot1.ceu.vcf";
    private static final String baseCommandPrefix = "-T IndelRealigner -noPG -R " + b36KGReference + " -I " + mainTestBam + " -targetIntervals " + mainTestIntervals + " -compress 0 -L 20:49,500-55,500 ";
    private static final String baseCommand = baseCommandPrefix + "-o %s ";
    private static final String base_md5 = "7574ab7d0b1ee5d44a0b3f85b6e944e6";
    private static final String base_md5_with_SW_or_VCF = "a918d69d26d3c87b29002ed31f428c48";

    @Test
    public void testDefaults() {

        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand,
                1,
                Arrays.asList(base_md5));
        executeTest("test realigner defaults", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommand + "-known " + knownIndels,
                1,
                Arrays.asList(base_md5_with_SW_or_VCF));
        executeTest("test realigner defaults with VCF", spec2);
    }

    @Test
    public void testKnownsOnly() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand + "--consensusDeterminationModel KNOWNS_ONLY -known " + knownIndels,
                1,
                Arrays.asList("36718f10d523dfb0fa2a709480f24bd4"));
        executeTest("realigner known indels only from VCF", spec1);
    }

    @Test
    public void testUseSW() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand + "--consensusDeterminationModel USE_SW -known " + knownIndels,
                1,
                Arrays.asList(base_md5_with_SW_or_VCF));
        executeTest("realigner use SW from VCF", spec1);
    }

    @Test
    public void testLods() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "-LOD 60", base_md5 );
        e.put( "-LOD 1 --consensusDeterminationModel USE_SW",  "9a75a0f7ad0442c78d0f8df260e733a4" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + entry.getKey(),
                    1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("realigner [%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testLongRun() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10,000,000-11,000,000 -targetIntervals " + validationDataLocation + "indelRealignerTest.NA12878.chrom1.intervals -compress 0 -o %s",
                1,
                Arrays.asList("e98f51d71f0a82141b36a7e9f94db237"));
        executeTest("realigner long run", spec);
    }

    @Test
    public void testNoTags() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand + "--noOriginalAlignmentTags --consensusDeterminationModel USE_SW",
                1,
                Arrays.asList("58ac675d0699eb236d469b8e84513d11"));
        executeTest("realigner no output tags", spec);
    }

    @Test
    public void testStats() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommandPrefix + "-stats %s -o /dev/null",
                1,
                Arrays.asList("7ed8d4eed635613fd031598a5c9ef5a3"));
        executeTest("realigner stats", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommandPrefix + "-LOD 60 -stats %s -o /dev/null",
                1,
                Arrays.asList("e8b02bfc5debec55fe936a38c59463cc"));
        executeTest("realigner stats", spec2);
    }

    @Test
    public void testMaxReadsInMemory() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "--maxReadsInMemory 10000", base_md5 );
        e.put( "--maxReadsInMemory 40000", base_md5 );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + entry.getKey(),
                    1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("realigner [%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testNWayOut() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommandPrefix + " -nWayOut .clean.bam ",
                1,
                Arrays.asList("d41d8cd98f00b204e9800998ecf8427e"));
        executeTest("test realigner nWayOut", spec1);
    }

}

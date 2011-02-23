package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class IndelRealignerIntegrationTest extends WalkerTest {

    private static final String mainTestBam = validationDataLocation + "indelRealignerTest.pilot1.ceu.bam";
    private static final String mainTestIntervals = validationDataLocation + "indelRealignerTest.pilot1.ceu.intervals";
    private static final String knownIndels = validationDataLocation + "indelRealignerTest.pilot1.ceu.vcf";
    private static final String baseCommandPrefix = "-T IndelRealigner -noPG -R " + b36KGReference + " -I " + mainTestBam + " -targetIntervals " + mainTestIntervals + " -compress 0 -L 20:49,500-55,500 --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe ";
    private static final String baseCommand = baseCommandPrefix + "-o %s ";

    @Test
    public void testDefaults() {
        String md5 = "20ff8b76d834a8aaca46405e8328d258";

        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand,
                1,
                Arrays.asList(md5));
        executeTest("test realigner defaults", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommand + "-B:indels,vcf " + knownIndels,
                1,
                Arrays.asList(md5));
        executeTest("test realigner defaults with VCF", spec2);

        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseCommand + "-D " + GATKDataLocation + "dbsnp_129_b36.rod",
                1,
                Arrays.asList(md5));
        executeTest("realigner defaults with dbsnp", spec3);

    }

    @Test
    public void testStats() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommandPrefix + "-stats %s -o /dev/null",
                1,
                Arrays.asList("ed5a207ddf5bdda4bb76899fb3eae35c"));
        executeTest("realigner stats", spec);

    }

    @Test
    public void testLods() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "-LOD 60", "20ff8b76d834a8aaca46405e8328d258" );
        e.put( "-LOD 1", "39862f48d9eaaf841ca4a0d2c05c4187" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + entry.getKey(),
                    1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("realigner [%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testKnownsOnly() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand + "-knownsOnly -B:indels,vcf " + knownIndels,
                1,
                Arrays.asList("1218d3c8fbd50581af5815938d6c0070"));
        executeTest("realigner known indels only from VCF", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommand + "-knownsOnly -D " + GATKDataLocation + "dbsnp_129_b36.rod",
                1,
                Arrays.asList("848740b201adc5c45bf82384c1f19d4d"));
        executeTest("realigner known indels only from dbsnp", spec2);
    }

    @Test
    public void testNoTags() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand + "--noOriginalAlignmentTags",
                1,
                Arrays.asList("00e009c97905f4fa89e3102261a1fd57"));
        executeTest("realigner no output tags", spec);
    }

    @Test
    public void testLongRun() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10,000,000-11,000,000 -targetIntervals " + validationDataLocation + "indelRealignerTest.NA12878.chrom1.intervals -compress 0 --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe -o %s",
                1,
                Arrays.asList(""));
        executeTest("realigner long run", spec);
    }
}

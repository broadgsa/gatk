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
    private static final String baseCommandPrefix = "-T IndelRealigner -noPG -R " + b36KGReference + " -I " + mainTestBam + " -targetIntervals " + mainTestIntervals + " -compress 0 -L 20:49,500-55,500 ";
    private static final String baseCommand = baseCommandPrefix + "-o %s ";
    private static final String base_md5 = "282070822dc5495eb20dad157d827133";

    @Test
    public void testDefaults() {

        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand,
                1,
                Arrays.asList(base_md5));
        executeTest("test realigner defaults", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommand + "-B:indels,vcf " + knownIndels,
                1,
                Arrays.asList(base_md5));
        executeTest("test realigner defaults with VCF", spec2);

        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseCommand + "-D " + GATKDataLocation + "dbsnp_129_b36.rod",
                1,
                Arrays.asList(base_md5));
        executeTest("realigner defaults with dbsnp", spec3);

    }

    @Test(dependsOnMethods = { "testNoTags" })
    public void testStats() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommandPrefix + "-stats %s -o /dev/null",
                1,
                Arrays.asList("ed5a207ddf5bdda4bb76899fb3eae35c"));
        executeTest("realigner stats", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommandPrefix + "-LOD 60 -stats %s -o /dev/null",
                1,
                Arrays.asList("ffab7d9ca19daa8a21e0b8f0072d39e9"));
        executeTest("realigner stats", spec2);        
    }

    @Test(dependsOnMethods = { "testKnownsOnly" })
    public void testLods() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "-LOD 60", base_md5 );
        e.put( "-LOD 1",  "c98d699d94f01bd0089f12646c764dfc" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + entry.getKey(),
                    1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("realigner [%s]", entry.getKey()), spec);
        }
    }

    @Test(dependsOnMethods = { "testDefaults" })
    public void testKnownsOnly() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseCommand + "-knownsOnly -B:indels,vcf " + knownIndels,
                1,
                Arrays.asList("36644c80f5e7b7c8679c0485ef681cd8"));
        executeTest("realigner known indels only from VCF", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseCommand + "-knownsOnly -D " + GATKDataLocation + "dbsnp_129_b36.rod",
                1,
                Arrays.asList("eab2cce434435da7dabb0926101c5586"));
        executeTest("realigner known indels only from dbsnp", spec2);
    }

    @Test(dependsOnMethods = { "testLongRun" })
    public void testNoTags() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand + "--noOriginalAlignmentTags",
                1,
                Arrays.asList("00ecb9df5afe3e9d61a75a2d019cb425"));
        executeTest("realigner no output tags", spec);
    }

    @Test(dependsOnMethods = { "testLods" })
    public void testLongRun() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10,000,000-11,000,000 -targetIntervals " + validationDataLocation + "indelRealignerTest.NA12878.chrom1.intervals -compress 0 -o %s",
                1,
                Arrays.asList("be859f9a98d738becee0526887cae42e"));
        executeTest("realigner long run", spec);
    }
}

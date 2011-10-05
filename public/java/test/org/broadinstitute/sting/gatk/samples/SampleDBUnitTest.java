package org.broadinstitute.sting.gatk.samples;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Sep 9, 2010
 * Time: 8:21:00 AM
 */
public class SampleDBUnitTest extends BaseTest {
    private static SampleDBBuilder builder;
    // all the test sample files are located here
    private File testPED = new File(testDir +  "ceutrio.ped");

    private static final Set<Sample> testPEDSamples = new HashSet<Sample>(Arrays.asList(
            new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED),
            new Sample("dad", "fam1", null, null,   Gender.MALE,   Affection.UNAFFECTED),
            new Sample("mom", "fam1", null, null,   Gender.FEMALE, Affection.AFFECTED)));

    private static final Set<Sample> testSAMSamples = new HashSet<Sample>(Arrays.asList(
            new Sample("kid", null, null, null, Gender.UNKNOWN,   Affection.UNKNOWN),
            new Sample("mom", null, null, null, Gender.UNKNOWN,   Affection.UNKNOWN),
            new Sample("dad", null, null, null, Gender.UNKNOWN,   Affection.UNKNOWN)));

    private static final String testPEDString =
            String.format("%s%n%s%n%s",
                    "fam1 kid dad mom 1 2",
                    "fam1 dad 0   0   1 1",
                    "fam1 mom 0   0   2 2");

    private static final String testPEDMultipleFamilies =
            String.format("%s%n%s%n%s%n%s%n%s",
                    "fam1 kid dad mom 1 2",
                    "fam1 dad 0   0   1 1",
                    "fam1 mom 0   0   2 2",
                    "fam3 s1  d1  m1  2 2",
                    "fam2 s2  d2  m2  2 2");

    private static final String testPEDStringInconsistentGender =
            "fam1 kid 0   0   2 2";

    private static final Set<Sample> testPEDSamplesAsSet =
            new HashSet<Sample>(testPEDSamples);


    @BeforeMethod
    public void before() {
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT);
    }

    @Test()
    public void loadPEDFile() {
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamplesAsSet, db.getSamples());
    }

    @Test()
    public void loadPEDString() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDString));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamplesAsSet, db.getSamples());
    }

    private static final void addSAMHeader() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        ArtificialSAMUtils.createEnumeratedReadGroups(header, Arrays.asList("1", "2", "3"),
                Arrays.asList("kid", "mom", "dad"));
        builder.addSamplesFromSAMHeader(header);
    }

    @Test()
    public void loadSAMHeader() {
        addSAMHeader();
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testSAMSamples, db.getSamples());
    }

    @Test()
    public void loadSAMHeaderPlusPED() {
        addSAMHeader();
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamples, db.getSamples());
    }

    @Test()
    public void loadDuplicateData() {
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamples, db.getSamples());
    }

    @Test(expectedExceptions = UserException.class)
    public void loadNonExistentFile() {
        builder.addSamplesFromPedigreeFiles(Arrays.asList(new File("non-existence-file.txt")));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testSAMSamples, db.getSamples());
    }

    @Test(expectedExceptions = UserException.class)
    public void loadInconsistentData() {
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringInconsistentGender));
        builder.getFinalSampleDB();
    }

    @Test(expectedExceptions = UserException.class)
    public void sampleInSAMHeaderNotInSamplesDB() {
        addSAMHeader();
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringInconsistentGender));
        builder.getFinalSampleDB();
    }

    @Test()
    public void getFamilyIDs() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getFamilyIDs(), new TreeSet<String>(Arrays.asList("fam1", "fam2", "fam3")));
    }

    @Test()
    public void getFamily() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getFamily("fam1"), testPEDSamplesAsSet);
    }

    @Test()
    public void loadFamilyIDs() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Map<String, Set<Sample>> families = db.getFamilies();
        Assert.assertEquals(families.size(), 3);
        Assert.assertEquals(families.keySet(), new TreeSet<String>(Arrays.asList("fam1", "fam2", "fam3")));

        for ( final String famID : families.keySet() ) {
            final Set<Sample> fam = families.get(famID);
            Assert.assertEquals(fam.size(), 3);
            for ( final Sample sample : fam ) {
                Assert.assertEquals(sample.getFamilyID(), famID);
            }
        }
    }
}

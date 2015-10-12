/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.samples;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
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
    private File testPED = new File(privateTestDir +  "testtrio.ped");

    private static final Set<Sample> testPEDSamples = new HashSet<Sample>(Arrays.asList(
            new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED),
            new Sample("dad", "fam1", null, null,   Gender.MALE,   Affection.UNAFFECTED),
            new Sample("mom", "fam1", null, null,   Gender.FEMALE, Affection.AFFECTED)));

    private static final Set<Sample> testPEDFamilyF2 = new HashSet<Sample>(Arrays.asList(
            new Sample("s2", "fam2", "d2", "m2", Gender.FEMALE, Affection.AFFECTED),
            new Sample("d2", "fam2", null, null, Gender.MALE, Affection.UNKNOWN),
            new Sample("m2", "fam2", null, null, Gender.FEMALE, Affection.UNKNOWN)
            ));

    private static final Set<Sample> testPEDFamilyF3 = new HashSet<Sample>(Arrays.asList(
            new Sample("s1", "fam3", "d1", "m1", Gender.FEMALE, Affection.AFFECTED),
            new Sample("d1", "fam3", null, null, Gender.MALE, Affection.UNKNOWN),
            new Sample("m1", "fam3", null, null, Gender.FEMALE, Affection.UNKNOWN)
            ));

    private static final Set<Sample> testSAMSamples = new HashSet<Sample>(Arrays.asList(
            new Sample("kid", null, null, null, Gender.UNKNOWN,   Affection.UNKNOWN),
            new Sample("mom", null, null, null, Gender.UNKNOWN,   Affection.UNKNOWN),
            new Sample("dad", null, null, null, Gender.UNKNOWN,   Affection.UNKNOWN)));

    private static final HashMap<String, Set<Sample>> testGetFamilies = new HashMap<String,Set<Sample>>();
    static {
        testGetFamilies.put("fam1", testPEDSamples);
        testGetFamilies.put("fam2", testPEDFamilyF2);
        testGetFamilies.put("fam3", testPEDFamilyF3);
    }

    private static final Set<Sample> testKidsWithParentsFamilies2 = new HashSet<Sample>(Arrays.asList(
            new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED),
            new Sample("kid3", "fam5", "dad2", "mom2", Gender.MALE,   Affection.AFFECTED),
            new Sample("kid2", "fam5", "dad2", "mom2", Gender.MALE,   Affection.AFFECTED)));

    private static final HashSet<String> testGetPartialFamiliesIds =   new HashSet<String>(Arrays.asList("kid","s1"));
    private static final HashMap<String, Set<Sample>> testGetPartialFamilies = new HashMap<String,Set<Sample>>();
    static {
        testGetPartialFamilies.put("fam1", new HashSet<Sample>(Arrays.asList(new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED))));
        testGetPartialFamilies.put("fam3", new HashSet<Sample>(Arrays.asList(new Sample("s1", "fam3", "d1", "m1", Gender.FEMALE, Affection.AFFECTED))));
    }

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

    private static final String testPEDMultipleFamilies2 =
            String.format("%s%n%s%n%s%n%s%n%s%n%s%n%s%n%s%n%s",
                    "fam1 kid dad mom 1 2",
                    "fam1 dad 0   0   1 1",
                    "fam1 mom 0   0   2 2",
                    "fam4 kid4 dad4 0 1 2",
                    "fam4 dad4 0   0   1 1",
                    "fam5 kid2 dad2 mom2 1 2",
                    "fam5 kid3 dad2 mom2 1 2",
                    "fam5 dad2 0   0   1 1",
                    "fam5 mom2 0   0   2 2");

    private static final String testPEDStringInconsistentGender =
            "fam1 kid 0   0   2 2";

    private static final String testPEDStringConsistent =
            "fam1 kid dad   mom   1 2";

    private static final Set<Sample> testPEDSamplesAsSet =
            new HashSet<Sample>(testPEDSamples);


    @BeforeMethod
    public void before() {
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT);
    }

    @Test()
    public void loadPEDFile() {
        final SampleDB db = builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED))
                                   .getFinalSampleDB();
        Assert.assertEquals(testPEDSamplesAsSet, db.getSamples());
    }

    @Test()
    public void loadPEDString() {
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDString))
                             .getFinalSampleDB();
        Assert.assertEquals(testPEDSamplesAsSet, db.getSamples());
    }

    private static final void addSAMHeader() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        ArtificialSAMUtils.createEnumeratedReadGroups(header, Arrays.asList("1", "2", "3"),
                Arrays.asList("kid", "mom", "dad"));
        builder.addSamplesFromSAMHeader(header);
    }

    @Test()
    public void loadSAMHeader() {
        addSAMHeader();
        final SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testSAMSamples, db.getSamples());
    }

    @Test()
    public void loadSAMHeaderPlusPED() {
        addSAMHeader();
        final SampleDB db = builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED))
                                   .getFinalSampleDB();
        Assert.assertEquals(testPEDSamples, db.getSamples());
    }

    @Test()
    public void loadDuplicateData() {
        final SampleDB db = builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED))
                                   .addSamplesFromPedigreeFiles(Arrays.asList(testPED))
                                   .getFinalSampleDB();
        Assert.assertEquals(testPEDSamples, db.getSamples());
    }

    @Test(expectedExceptions = UserException.class)
    public void loadNonExistentFile() {
        final SampleDB db = builder.addSamplesFromPedigreeFiles(Arrays.asList(new File("non-existence-file.txt")))
                           .getFinalSampleDB();
        Assert.assertEquals(testSAMSamples, db.getSamples());
    }

    @Test(expectedExceptions = UserException.class)
    public void loadInconsistentData() {
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT)
                      .addSamplesFromPedigreeFiles(Arrays.asList(testPED))
                      .addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringInconsistentGender));
        builder.getFinalSampleDB();
    }

    @Test
    public void loadConsistentData() {
        // build a temporary DB and get the resulting sample to use for test result comparison
        final Sample baseKidSample = new SampleDBBuilder(PedigreeValidationType.STRICT)
                                        .addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringConsistent))
                                        .getFinalSampleDB()
                                        .getSample("kid");

        // build a sample DB and then merge in the consistent test string
        final SampleDB finalDB = new SampleDBBuilder(PedigreeValidationType.STRICT)
                                     .addSamplesFromPedigreeFiles(Arrays.asList(testPED))
                                     .addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringConsistent))
                                     .getFinalSampleDB();

        Assert.assertEquals(finalDB.getSamples().size(), 3);
        Assert.assertTrue(finalDB.getSample("kid").equals(baseKidSample));
    }

    @Test(expectedExceptions = UserException.class)
    public void sampleInSAMHeaderNotInSamplesDB() {
        addSAMHeader();
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringInconsistentGender))
               .getFinalSampleDB();
    }

    @Test()
    public void getFamilyIDs() {
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies))
                                   .getFinalSampleDB();
        Assert.assertEquals(db.getFamilyIDs(), new TreeSet<String>(Arrays.asList("fam1", "fam2", "fam3")));
    }

    @Test()
    public void getFamily() {
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies))
                                   .getFinalSampleDB();
        Assert.assertEquals(db.getFamily("fam1"), testPEDSamplesAsSet);
    }

    @Test()
    public void getFamilies(){
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies))
                                   .getFinalSampleDB();
        Assert.assertEquals(db.getFamilies(),testGetFamilies);
        Assert.assertEquals(db.getFamilies(null),testGetFamilies);
        Assert.assertEquals(db.getFamilies(testGetPartialFamiliesIds),testGetPartialFamilies);
    }

    @Test()
    public void testGetChildrenWithParents()
    {
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies2))
                                   .getFinalSampleDB();
        Assert.assertEquals(db.getChildrenWithParents(), testKidsWithParentsFamilies2);
        Assert.assertEquals(db.getChildrenWithParents(false), testKidsWithParentsFamilies2);
        Assert.assertEquals(db.getChildrenWithParents(true), new HashSet<Sample>(Arrays.asList(new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED))));
    }

    @Test()
    public void testGetFounderIds(){
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies2))
                                   .getFinalSampleDB();
        Assert.assertEquals(db.getFounderIds(), new HashSet<String>(Arrays.asList("dad","mom","dad2","mom2","dad4")));
    }

    @Test()
    public void loadFamilyIDs() {
        final SampleDB db = builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies))
                                   .getFinalSampleDB();
        final Map<String, Set<Sample>> families = db.getFamilies();
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

/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.samples;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.StringReader;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

/**
 * UnitTest for PedReader
 *
 * @author Mark DePristo
 * @since 2011
 */
public class PedReaderUnitTest extends BaseTest {
    private static Logger logger = Logger.getLogger(PedReaderUnitTest.class);

    private class PedReaderTest extends TestDataProvider {
        public String fileContents;
        public List<Sample> expectedSamples;

        private PedReaderTest(final String name, final List<Sample> expectedSamples, final String fileContents) {
            super(PedReaderTest.class, name);
            this.fileContents = fileContents;
            this.expectedSamples = expectedSamples;
        }
    }

//     Family ID
//     Individual ID
//     Paternal ID
//     Maternal ID
//     Sex (1=male; 2=female; other=unknown)
//     Phenotype
//
//     -9 missing
//     0 missing
//     1 unaffected
//     2 affected

    @DataProvider(name = "readerTest")
    public Object[][] createPEDFiles() {
        new PedReaderTest("singleRecordMale",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.MALE, Affection.UNAFFECTED)),
                "fam1 kid 0 0 1 1");

        new PedReaderTest("singleRecordFemale",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.FEMALE, Affection.UNAFFECTED)),
                "fam1 kid 0 0 2 1");

        new PedReaderTest("singleRecordMissingGender",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.UNKNOWN, Affection.UNKNOWN)),
                "fam1 kid 0 0 0 0");

        // Affection
        new PedReaderTest("singleRecordAffected",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.MALE, Affection.AFFECTED)),
                "fam1 kid 0 0 1 2");

        new PedReaderTest("singleRecordUnaffected",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.MALE, Affection.UNAFFECTED)),
                "fam1 kid 0 0 1 1");

        new PedReaderTest("singleRecordMissingAffection-9",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.MALE, Affection.UNKNOWN)),
                "fam1 kid 0 0 1 -9");

        new PedReaderTest("singleRecordMissingAffection0",
                Arrays.asList(new Sample("kid", "fam1", null, null, Gender.MALE, Affection.UNKNOWN)),
                "fam1 kid 0 0 1 0");

        new PedReaderTest("multipleUnrelated",
                Arrays.asList(
                    new Sample("s1", "fam1", null, null, Gender.MALE,   Affection.AFFECTED),
                    new Sample("s2", "fam2", null, null, Gender.FEMALE, Affection.UNAFFECTED)),
                String.format("%s\n%s",
                        "fam1 s1 0 0 1 1",
                        "fam2 s2 0 0 2 2"));

        new PedReaderTest("explicitTrio",
                Arrays.asList(
                    new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED),
                    new Sample("dad", "fam1", null, null,   Gender.MALE,   Affection.UNAFFECTED),
                    new Sample("mom", "fam1", null, null,   Gender.FEMALE, Affection.AFFECTED)),
                String.format("%s\n%s\n%s",
                        "fam1 kid dad mom 1 2",
                        "fam1 dad 0   0   1 1",
                        "fam1 mom 0   0   2 2"));

        new PedReaderTest("implicitTrio",
                Arrays.asList(
                    new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED),
                    new Sample("dad", "fam1", null, null,   Gender.MALE,   Affection.UNKNOWN),
                    new Sample("mom", "fam1", null, null,   Gender.FEMALE, Affection.UNKNOWN)),
                "fam1 kid dad mom 1 1");

        new PedReaderTest("partialTrio",
                Arrays.asList(
                    new Sample("kid", "fam1", "dad", "mom", Gender.MALE,   Affection.AFFECTED),
                    new Sample("dad", "fam1", null, null,   Gender.MALE,   Affection.UNAFFECTED),
                    new Sample("mom", "fam1", null, null,   Gender.FEMALE, Affection.UNKNOWN)),
                String.format("%s\n%s",
                        "fam1 kid dad mom 1 2",
                        "fam1 dad 0   0   1 1"));

        new PedReaderTest("bigPedigree",
                Arrays.asList(
                    new Sample("kid", "fam1", "dad",       "mom",      Gender.MALE,   Affection.AFFECTED),
                    new Sample("dad", "fam1", "granddad1", "grandma1", Gender.MALE,   Affection.UNAFFECTED),
                    new Sample("granddad1", "fam1", null, null,        Gender.MALE,   Affection.UNKNOWN),
                    new Sample("grandma1",  "fam1", null, null,        Gender.FEMALE,   Affection.UNKNOWN),
                    new Sample("mom", "fam1", "granddad2", "grandma2", Gender.FEMALE, Affection.AFFECTED),
                    new Sample("granddad2", "fam1", null, null,        Gender.MALE,   Affection.UNKNOWN),
                    new Sample("grandma2",  "fam1", null, null,        Gender.FEMALE,   Affection.UNKNOWN)),
                String.format("%s\n%s\n%s",
                        "fam1 kid dad       mom      1 2",
                        "fam1 dad granddad1 grandma1 1 1",
                        "fam1 mom granddad2 grandma2 2 2"));

        // Quantitative trait
        new PedReaderTest("QuantitativeTrait",
                Arrays.asList(
                    new Sample("s1", "fam1", null, null, Gender.MALE,   Affection.QUANTITATIVE, 1.0),
                    new Sample("s2", "fam2", null, null, Gender.FEMALE, Affection.QUANTITATIVE, 10.0)),
                String.format("%s\n%s",
                        "fam1 s1 0 0 1 1",
                        "fam2 s2 0 0 2 10.0"));

        new PedReaderTest("QuantitativeTraitWithMissing",
                Arrays.asList(
                    new Sample("s1", "fam1", null, null, Gender.MALE,   Affection.UNKNOWN, Sample.UNSET_QT),
                    new Sample("s2", "fam2", null, null, Gender.FEMALE, Affection.QUANTITATIVE, 10.0)),
                String.format("%s\n%s",
                        "fam1 s1 0 0 1 -9",
                        "fam2 s2 0 0 2 10.0"));

        new PedReaderTest("QuantitativeTraitOnlyInts",
                Arrays.asList(
                    new Sample("s1", "fam1", null, null, Gender.MALE,   Affection.QUANTITATIVE, 1.0),
                    new Sample("s2", "fam2", null, null, Gender.FEMALE, Affection.QUANTITATIVE, 10.0)),
                String.format("%s\n%s",
                        "fam1 s1 0 0 1 1",
                        "fam2 s2 0 0 2 10"));

        return PedReaderTest.getTests(PedReaderTest.class);
    }

    private static final void runTest(PedReaderTest test, String myFileContents, EnumSet<PedReader.MissingPedFields> missing) {
        logger.warn("Test " + test);
        PedReader reader = new PedReader();
        SampleDataSource sampleDB = new SampleDataSource();
        List<Sample> readSamples = reader.parse(new StringReader(myFileContents), missing, sampleDB);
        Assert.assertEquals(test.expectedSamples, readSamples, "Parsed incorrect number of samples");
    }

    @Test(enabled = true, dataProvider = "readerTest")
    public void testPedReader(PedReaderTest test) {
        runTest(test, test.fileContents, EnumSet.noneOf(PedReader.MissingPedFields.class));
    }

    @Test(enabled = true, dataProvider = "readerTest", dependsOnMethods = "testPedReader")
    public void testPedReaderWithComments(PedReaderTest test) {
        runTest(test, "#comment\n" + test.fileContents, EnumSet.noneOf(PedReader.MissingPedFields.class));
    }

    @Test(enabled = true, dataProvider = "readerTest", dependsOnMethods = "testPedReader")
    public void testPedReaderWithMissing(PedReaderTest test) {
        // todo -- test MISSING by splicing strings
        //runTest(test, "#comment\n" + test.fileContents, EnumSet.noneOf(PedReader.MissingPedFields.class));
    }

}
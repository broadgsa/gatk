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

package org.broadinstitute.gatk.utils.downsampling;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;


/**
 * Basic unit test for AlleleBiasedDownsamplingUtils
 */
public class AlleleBiasedDownsamplingUtilsUnitTest extends BaseTest {


    @Test
    public void testSmartDownsampling() {

        final int[] idealHetAlleleCounts = new int[]{0, 50, 0, 50};
        final int[] idealHomAlleleCounts = new int[]{0, 100, 0, 0};

        // no contamination, no removal
        testOneCase(0, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);

        // hom sample, het contaminant, different alleles
        testOneCase(5, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 5, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 0, 5, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);

        // hom sample, hom contaminant, different alleles
        testOneCase(10, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 10, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 0, 10, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);

        // het sample, het contaminant, different alleles
        testOneCase(5, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 5, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);

        // het sample, hom contaminant, different alleles
        testOneCase(10, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 10, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);

        // hom sample, het contaminant, overlapping alleles
        final int[] enhancedHomAlleleCounts = new int[]{0, 105, 0, 0};
        testOneCase(5, 5, 0, 0, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts);
        testOneCase(0, 5, 5, 0, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts);
        testOneCase(0, 5, 0, 5, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts);

        // hom sample, hom contaminant, overlapping alleles
        testOneCase(0, 10, 0, 0, 0.1, 100, idealHomAlleleCounts, new int[]{0, 110, 0, 0});

        // het sample, het contaminant, overlapping alleles
        testOneCase(5, 5, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 5, 5, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 5, 0, 5, 0.1, 100, idealHetAlleleCounts, new int[]{0, 55, 0, 55});
        testOneCase(5, 0, 0, 5, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 5, 5, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);

        // het sample, hom contaminant, overlapping alleles
        testOneCase(0, 10, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 0, 10, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
    }

    private static void testOneCase(final int addA, final int addC, final int addG, final int addT, final double contaminationFraction,
                                    final int pileupSize, final int[] initialCounts, final int[] targetCounts) {

        final int[] actualCounts = initialCounts.clone();
        actualCounts[0] += addA;
        actualCounts[1] += addC;
        actualCounts[2] += addG;
        actualCounts[3] += addT;

        final int[] results = AlleleBiasedDownsamplingUtils.runSmartDownsampling(actualCounts, (int) (pileupSize * contaminationFraction));
        Assert.assertTrue(countsAreEqual(results, targetCounts));
    }

    private static boolean countsAreEqual(final int[] counts1, final int[] counts2) {
        for ( int i = 0; i < 4; i++ ) {
            if ( counts1[i] != counts2[i] )
                return false;
        }
        return true;
    }

    @DataProvider(name = "BiasedDownsamplingTest")
    public Object[][] makeBiasedDownsamplingTest() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

        for ( final int originalCount : Arrays.asList(1, 2, 10, 1000) ) {
            for ( final int toRemove : Arrays.asList(0, 1, 2, 10, 1000) ) {
                if ( toRemove <= originalCount )
                    tests.add(new Object[]{header, originalCount, toRemove});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BiasedDownsamplingTest")
    public void testBiasedDownsampling(final SAMFileHeader header, final int originalCount, final int toRemove) {

        final LinkedList<PileupElement> elements = new LinkedList<>();
        for ( int i = 0; i < originalCount; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, 1, 1);
            elements.add(new PileupElement(read, 0, new CigarElement(1, CigarOperator.M), 0, 0));
        }

        final List<PileupElement> result = AlleleBiasedDownsamplingUtils.downsampleElements(elements, originalCount, toRemove);

        Assert.assertEquals(result.size(), toRemove);
    }

    @Test
    public void testLoadContaminationFileDetails(){
        Logger logger=org.apache.log4j.Logger.getRootLogger();

        final String ArtificalBAMLocation = privateTestDir + "ArtificallyContaminatedBams/";
        final File ContamFile1=new File(ArtificalBAMLocation+"contamination.case.1.txt");

        Map<String,Double> Contam1=new HashMap<String,Double>();
        Set<String> Samples1=new HashSet<String>();

        Contam1.put("NA11918",0.15);
        Samples1.addAll(Contam1.keySet());
        testLoadFile(ContamFile1,Samples1,Contam1,logger);

        Contam1.put("NA12842",0.13);
        Samples1.addAll(Contam1.keySet());
        testLoadFile(ContamFile1,Samples1,Contam1,logger);

        Samples1.add("DUMMY");
        testLoadFile(ContamFile1,Samples1,Contam1,logger);
   }

    private static void testLoadFile(final File file, final Set<String> Samples, final Map<String,Double> map, Logger logger){
        Map<String,Double> loadedMap = AlleleBiasedDownsamplingUtils.loadContaminationFile(file,0.0,Samples,logger);
        Assert.assertTrue(loadedMap.equals(map));
    }

    @DataProvider(name = "goodContaminationFiles")
    public Integer[][] goodContaminationFiles() {
        return new Integer[][]{
                {1, 2},
                {2, 3},
                {3, 2},
                {4, 2},
                {5, 3},
                {6, 2},
                {7, 2},
                {8, 2}
        };
    }

    @Test(dataProvider = "goodContaminationFiles")
    public void testLoadContaminationFile(final Integer ArtificalBAMnumber, final Integer numberOfSamples) {
        final String ArtificialBAM = String.format("ArtificallyContaminatedBams/contamination.case.%d.txt", ArtificalBAMnumber);
        Logger logger = org.apache.log4j.Logger.getRootLogger();

        File ContamFile = new File(privateTestDir, ArtificialBAM);
        Assert.assertTrue(AlleleBiasedDownsamplingUtils.loadContaminationFile(ContamFile, 0.0, null, logger).size() == numberOfSamples);

    }


    @DataProvider(name = "badContaminationFiles")
    public Integer[][] badContaminationFiles() {
        return new Integer[][]{{1}, {2}, {3}, {4}, {5}};
    }

    @Test(dataProvider = "badContaminationFiles", expectedExceptions = UserException.MalformedFile.class)
    public void testLoadBrokenContaminationFile(final int i) {
        Logger logger = org.apache.log4j.Logger.getRootLogger();
        final String ArtificalBAMLocation = privateTestDir + "ArtificallyContaminatedBams/";

        File ContaminationFile = new File(ArtificalBAMLocation + String.format("contamination.case.broken.%d.txt", i));
        AlleleBiasedDownsamplingUtils.loadContaminationFile(ContaminationFile, 0.0, null, logger);

    }


}

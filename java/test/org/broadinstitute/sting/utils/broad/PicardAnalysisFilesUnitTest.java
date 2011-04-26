package org.broadinstitute.sting.utils.broad;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;

import static org.broadinstitute.sting.utils.broad.PicardAggregationUtilsUnitTest.*;

public class PicardAnalysisFilesUnitTest extends BaseTest {
    @Test
    public void testParseLatest() throws Exception {
        PicardAnalysisFiles files = new PicardAnalysisFiles(PROJECT, SAMPLE);
        Assert.assertNotNull(files.getPath());
        files = new PicardAnalysisFiles(PROJECT, SAMPLE, PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE));
        Assert.assertNotNull(files.getPath());
    }

    @Test
    public void testParseValid() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file.txt");
        Assert.assertEquals(file.getReferenceSequence(), "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
        Assert.assertEquals(file.getTargetIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list");
        Assert.assertEquals(file.getBaitIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list");
    }

    @Test
    public void testParseValidWithComments() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file_with_comments.txt");
        Assert.assertEquals(file.getReferenceSequence(), "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
        Assert.assertEquals(file.getTargetIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list");
        Assert.assertEquals(file.getBaitIntervals(), "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list");
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testParseBadPath() throws Exception {
        new PicardAnalysisFiles(BaseTest.validationDataLocation + "non_existent_picard_analysis_file.txt");
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testParseMissingLatest() throws Exception {
        new PicardAnalysisFiles(MISSING_PROJECT, MISSING_SAMPLE);
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testParseMissingVersion() throws Exception {
        new PicardAnalysisFiles(PROJECT, SAMPLE, PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE) + 2);
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testParseMultipleReferences() throws Exception {
        PicardAnalysisFiles file = new PicardAnalysisFiles(BaseTest.validationDataLocation + "picard_analysis_file_with_different_references.txt");
        file.getReferenceSequence();
    }
}

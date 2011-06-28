package org.broadinstitute.sting.utils.broad;

import org.testng.Assert;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.datasources.pipeline.Pipeline;
import org.broadinstitute.sting.datasources.pipeline.PipelineSample;
import org.broadinstitute.sting.utils.yaml.YamlUtils;
import org.testng.annotations.Test;
import static org.broadinstitute.sting.utils.broad.PicardAggregationUtilsUnitTest.*;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

public class PicardPipelineUnitTest {
    @Test
    public void testParseTsv() throws IOException {
        File tsv = writeTsv(PROJECT, SAMPLE);
        Pipeline pipeline = PicardPipeline.parse(tsv);
        validatePipeline(pipeline, FilenameUtils.getBaseName(tsv.getPath()));
    }

    @Test
    public void testParseTsvWithPicardComments() throws Exception {
        File tsv = writeTsv("C460", "HG01359");
        PicardPipeline.parse(tsv);
    }

    @Test
    public void testParseYaml() throws IOException {
        File yaml = writeYaml("project_name", PROJECT, SAMPLE);
        Pipeline pipeline = PicardPipeline.parse(yaml);
        validatePipeline(pipeline, "project_name");
    }

    private void validatePipeline(Pipeline pipeline, String name) {
        Assert.assertEquals(pipeline.getProject().getName(), name);
        Assert.assertTrue(pipeline.getProject().getReferenceFile().exists(), "reference not found");
        Assert.assertTrue(pipeline.getProject().getIntervalList().exists(), "intervals not found");
        Assert.assertTrue(pipeline.getProject().getRefseqTable().exists(), "refseq not found");
        Assert.assertTrue(pipeline.getProject().getGenotypeDbsnp().exists(), "genotype dbsnp not found");
        Assert.assertTrue(pipeline.getProject().getEvalDbsnp().exists(), "eval dbsnp not found");
        Assert.assertEquals(pipeline.getSamples().size(), 1);
        for (PipelineSample sample: pipeline.getSamples()) {
            Assert.assertEquals(sample.getId(), PROJECT + "_" + SAMPLE);
            Assert.assertTrue(sample.getBamFiles().get(PicardPipeline.PICARD_BAM_TYPE).exists(), "bam not found");
            Assert.assertEquals(sample.getTags().get(PicardPipeline.PROJECT_TAG), PROJECT);
            Assert.assertEquals(sample.getTags().get(PicardPipeline.SAMPLE_TAG), SAMPLE);
        }
    }

    private File writeTsv(String project, String sample) throws IOException {
        File tsv = BaseTest.createTempFile("pipeline", ".tsv");
        FileUtils.writeLines(tsv, Collections.singletonList(project + "\t" + sample));
        return tsv;
    }

    private File writeYaml(String projectName, String project, String sample) throws IOException {
        File yaml = BaseTest.createTempFile("pipeline", ".yaml");
        PipelineSample pipelineSample = new PipelineSample();
        pipelineSample.getTags().put(PicardPipeline.PROJECT_TAG, project);
        pipelineSample.getTags().put(PicardPipeline.SAMPLE_TAG, sample);
        Pipeline pipeline = new Pipeline();
        pipeline.getProject().setName(projectName);
        pipeline.getSamples().add(pipelineSample);
        YamlUtils.dump(pipeline, yaml);
        return yaml;
    }
}

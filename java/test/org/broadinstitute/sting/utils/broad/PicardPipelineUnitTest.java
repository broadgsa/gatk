package org.broadinstitute.sting.utils.broad;

import junit.framework.Assert;
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
    public void testParseYaml() throws IOException {
        File yaml = writeYaml("project_name", PROJECT, SAMPLE);
        Pipeline pipeline = PicardPipeline.parse(yaml);
        validatePipeline(pipeline, "project_name");
    }

    private void validatePipeline(Pipeline pipeline, String name) {
        Assert.assertEquals(pipeline.getProject().getName(), name);
        Assert.assertTrue("reference not found", pipeline.getProject().getReferenceFile().exists());
        Assert.assertTrue("intervals not found", pipeline.getProject().getIntervalList().exists());
        Assert.assertTrue("refseq not found", pipeline.getProject().getRefseqTable().exists());
        Assert.assertTrue("genotype dbsnp not found", pipeline.getProject().getGenotypeDbsnp().exists());
        Assert.assertTrue("eval dbsnp not found", pipeline.getProject().getEvalDbsnp().exists());
        Assert.assertEquals(pipeline.getSamples().size(), 1);
        for (PipelineSample sample: pipeline.getSamples()) {
            Assert.assertEquals(sample.getId(), PROJECT + "_" + SAMPLE);
            Assert.assertTrue("bam not found", sample.getBamFiles().get(PicardPipeline.PICARD_BAM_TYPE).exists());
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

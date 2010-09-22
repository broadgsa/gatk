/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.datasources.pipeline;

import org.broadinstitute.sting.utils.yaml.YamlUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.Map;

public class PipelineUnitTest {
    @Test
    public void testDumpAndLoad() throws Exception {
        Pipeline pipeline = new Pipeline();

        pipeline.getProject().setName("PRJ_NAME");
        pipeline.getProject().setReferenceFile(new File("my.fasta"));
        pipeline.getProject().setDbsnpFile(new File("my.dbsnp"));
        pipeline.getProject().getTags().put("testProjectTag", "project value here");

        PipelineSample sample = new PipelineSample();
        sample.setId("SMP_ID");
        sample.getBamFiles().put("recalibrated", new File("recalibrated.bam"));
        sample.getBamFiles().put("cleaned", new File("/absolute/path/to/cleaned.bam"));
        sample.getTags().put("testSampleTag", "sample value here");

        pipeline.getSamples().add(sample);

        File file = File.createTempFile("testDumpAndLoad", ".yaml");
        YamlUtils.dump(pipeline, file);
        Pipeline pipelineLoad = YamlUtils.load(Pipeline.class, file);

        Assert.assertEquals(pipeline.getProject().getName(), pipelineLoad.getProject().getName());
        Assert.assertEquals(pipeline.getProject().getReferenceFile(), pipelineLoad.getProject().getReferenceFile());
        Assert.assertEquals(pipeline.getProject().getIntervalList(), pipelineLoad.getProject().getIntervalList());
        Assert.assertEquals(pipeline.getProject().getDbsnpFile(), pipelineLoad.getProject().getDbsnpFile());

        Assert.assertEquals(pipeline.getProject().getTags().size(), pipelineLoad.getProject().getTags().size());
        for (Map.Entry<String, String> entry : pipeline.getProject().getTags().entrySet())
            Assert.assertEquals(entry.getValue(), pipeline.getProject().getTags().get(entry.getKey()));

        Assert.assertEquals(pipeline.getSamples().size(), pipelineLoad.getSamples().size());
        for (int i = 0; i < pipeline.getSamples().size(); i++) {
            PipelineSample pipelineSample = pipeline.getSamples().get(i);
            PipelineSample pipelineLoadSample = pipelineLoad.getSamples().get(i);

            Assert.assertEquals(pipelineSample.getId(), pipelineLoadSample.getId());

            Assert.assertEquals(pipelineSample.getBamFiles().size(), pipelineLoadSample.getBamFiles().size());
            for (Map.Entry<String, File> entry : pipelineSample.getBamFiles().entrySet())
                Assert.assertEquals(entry.getValue(), pipelineSample.getBamFiles().get(entry.getKey()));

            Assert.assertEquals(pipelineSample.getTags().size(), pipelineLoadSample.getTags().size());
            for (Map.Entry<String, String> entry : pipelineSample.getTags().entrySet())
                Assert.assertEquals(entry.getValue(), pipelineSample.getTags().get(entry.getKey()));
        }
    }
}

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

package org.broadinstitute.sting.utils.broad;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.NullArgumentException;
import org.broadinstitute.sting.datasources.pipeline.Pipeline;
import org.broadinstitute.sting.datasources.pipeline.PipelineProject;
import org.broadinstitute.sting.datasources.pipeline.PipelineSample;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.yaml.YamlUtils;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * Automatically gets the latest version using PicardAggregationUtils.
 */
public class PicardPipeline {

    protected static final String PROJECT_TAG = "SQUIDProject";
    protected static final String SAMPLE_TAG = "CollaboratorID";
    protected static final String PICARD_BAM_TYPE = "cleaned";

    private PicardPipeline() {}

    /**
     * Creates a new PicardPipeline
     * @param path Path to a tsv with project [tab] sample on each line or a pipeline yaml.
     * @return a new Picard
     * @throws FileNotFoundException when unable to find the file or any supporting files.
     */
    public static Pipeline parse(File path) throws FileNotFoundException {
        if (path == null)
            throw new NullArgumentException("path");

        Pipeline pipeline;
        if (path.getName().endsWith(".tsv")) {
            pipeline = new Pipeline();
            pipeline.getProject().setName(FilenameUtils.getBaseName(path.getPath()));
            for (String line: new XReadLines(path)) {
                String[] projectSample = line.split("\t");
                addSample(pipeline, projectSample[0], projectSample[1]);
            }
        } else if (path.getName().endsWith(".yaml")) {
            pipeline = YamlUtils.load(Pipeline.class, path);
        } else {
            throw new UserException.BadInput("Path does not end with .tsv or .yaml: " + path.getPath());
        }

        update(pipeline);
        return pipeline;
    }

    private static void update(Pipeline pipeline) throws FileNotFoundException {
        for (PipelineSample sample: pipeline.getSamples())
            updateSample(pipeline.getProject(), sample);
    }

    private static void addSample(Pipeline pipeline, String project, String sample) {
        PipelineSample pipelineSample = new PipelineSample();
        pipelineSample.getTags().put(PROJECT_TAG, project);
        pipelineSample.getTags().put(SAMPLE_TAG, sample);
        pipeline.getSamples().add(pipelineSample);
    }

    private static void updateSample(PipelineProject pipelineProject, PipelineSample pipelineSample) throws FileNotFoundException {
        if (!pipelineSample.getTags().containsKey(PROJECT_TAG) && !pipelineSample.getTags().containsKey(SAMPLE_TAG))
            return;
        String project = pipelineSample.getTags().get(PROJECT_TAG);
        String sample = pipelineSample.getTags().get(SAMPLE_TAG);
        int version = PicardAggregationUtils.getLatestVersion(project, sample);
        if (version <= 0)
            throw new UserException.BadInput("Project sample not found: " + project + "/" + sample);
        String bam = PicardAggregationUtils.getSampleBam(project, sample, version);
        if (pipelineSample.getId() == null)
            pipelineSample.setId(project + "_" + sample);
        pipelineSample.getBamFiles().put(PICARD_BAM_TYPE, new File(bam));

        PicardAnalysisFiles analysis = new PicardAnalysisFiles(project, sample, version);
        if (pipelineProject.getReferenceFile() == null) {
            String referenceSequence = analysis.getReferenceSequence();
            ReferenceData referenceData = ReferenceData.getByReference(referenceSequence);
            pipelineProject.setReferenceFile(new File(referenceData.getReference()));
            pipelineProject.setRefseqTable(new File(referenceData.getRefseq()));
            if (analysis.getTargetIntervals() != null)
                pipelineProject.setIntervalList(new File(analysis.getTargetIntervals()));
            pipelineProject.setEvalDbsnp(new File(referenceData.getDbsnp(129)));
            if (referenceData.getDbsnpVersions().contains(132)) {
                pipelineProject.setGenotypeDbsnp(new File(referenceData.getDbsnp(132)));
            } else {
                pipelineProject.setGenotypeDbsnp(new File(referenceData.getDbsnp(129)));
            }
        } else {
            String referenceSequence = analysis.getReferenceSequence();
            if (!pipelineProject.getReferenceFile().getPath().equals(referenceSequence))
                throw new UserException.BadInput("Samples sequenced with different references");
        }
    }
}

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

package org.broadinstitute.sting.pipeline;

import java.util.ArrayList;
import java.util.List;

/**
 * Java bean for storing a list of samples for a pipeline.
 *
 * NOTE: This class is used in a very similar way to the classes in
 * org.broadinstitute.sting.gatk.datasources.sample.
 *
 * Both store / load sample information from the file system as YAML.
 *
 * This package will likely be refactored to share common functionality
 * with the other at a future date as requirements coalesce.
 *
 * - kshakir September 22, 2010
 */
public class Pipeline {
    private PipelineProject project = new PipelineProject();
    private List<PipelineSample> samples = new ArrayList<PipelineSample>();

    public PipelineProject getProject() {
        return project;
    }

    public void setProject(PipelineProject project) {
        this.project = project;
    }

    public List<PipelineSample> getSamples() {
        return samples;
    }

    public void setSamples(List<PipelineSample> samples) {
        this.samples = samples;
    }
}

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

package org.broadinstitute.sting.queue.pipeline.playground

import org.testng.annotations.{DataProvider, Test}
import org.broadinstitute.sting.datasources.pipeline.{PipelineSample, Pipeline}
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.queue.pipeline._
import org.broadinstitute.sting.BaseTest

class HybridSelectionPipelineTest {
  def datasets = List(k1gChr20Dataset)

  val k1gChr20Dataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx_chr20", BaseTest.hg19Chr20Intervals)

    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.called.all.all.all", "nCalledLoci", 1392)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.called.all.known.all", "nCalledLoci", 1143)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.called.all.novel.all", "nCalledLoci", 249)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.called.all.all.all", "tiTvRatio", 3.6250)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.called.all.known.all", "tiTvRatio", 3.7190)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.called.all.novel.all", "tiTvRatio", 3.2037)

    dataset
  }

  def newK1gDataset(projectName: String, intervals: String) = {
    val project = PipelineTest.createHg19Project(projectName, intervals)
    var samples = List.empty[PipelineSample]
    for (k1gBam <- PipelineTest.k1gBams)
      samples :+= PipelineTest.createK1gSample(projectName, k1gBam)
    new PipelineDataset(PipelineTest.createPipeline(project, samples))
  }

  @DataProvider(name="datasets")//, parallel=true)
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets")
  def testHybridSelectionPipeline(dataset: PipelineDataset) {
    val projectName = dataset.pipeline.getProject.getName
    val testName = "HybridSelectionPipeline-" + projectName
    val yamlFile = writeYaml(testName, dataset.pipeline)

    // Run the pipeline with the expected inputs.
    val pipelineCommand =
      "-retry 1 -S scala/qscript/playground/HybridSelectionPipeline.scala -Y %s"
        .format(yamlFile)

    val pipelineSpec = new PipelineTestSpec
    pipelineSpec.name = testName
    pipelineSpec.args = pipelineCommand
    pipelineSpec.jobQueue = dataset.jobQueue

    pipelineSpec.evalSpec = new PipelineTestEvalSpec
    pipelineSpec.evalSpec.evalReport = projectName + ".eval"
    pipelineSpec.evalSpec.validations = dataset.validations

    PipelineTest.executeTest(pipelineSpec)
  }

  private def writeYaml(testName: String, pipeline: Pipeline) = {
    val yamlFile = BaseTest.createTempFile(pipeline.getProject.getName + ".", ".yaml")
    YamlUtils.dump(pipeline, yamlFile)
    yamlFile
  }

  class PipelineDataset(var pipeline: Pipeline = null,
                        var validations: List[PipelineValidation[_]] = Nil,
                        var jobQueue: String = null) {
    override def toString = pipeline.getProject.getName
  }
}

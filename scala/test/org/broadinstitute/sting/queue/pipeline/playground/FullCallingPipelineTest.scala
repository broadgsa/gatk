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
import collection.JavaConversions._
import java.io.File
import org.broadinstitute.sting.datasources.pipeline.{PipelineSample, Pipeline}
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.queue.pipeline._

class FullCallingPipelineTest {
  def datasets = List(k1gChr20Dataset)//, k1gExomeDataset)

  val k1gChr20Dataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx_chr20", true)

    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.all", "nCalledLoci", 1398)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.known", "nCalledLoci", 1143)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.novel", "nCalledLoci", 255)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.all", "tiTvRatio", 3.6250)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.known", "tiTvRatio", 3.7190)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.novel", "tiTvRatio", 3.2037)

    dataset
  }

  val k1gExomeDataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx", false)

    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.all", "nCalledLoci", 52668)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.known", "nCalledLoci", 41248)
    dataset.validations :+= new IntegerValidation("CountVariants", "dbsnp.eval.novel", "nCalledLoci", 11420)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.all", "tiTvRatio", 3.271)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.known", "tiTvRatio", 3.3299)
    dataset.validations :+= new DoubleValidation("TiTvVariantEvaluator", "dbsnp.eval.novel", "tiTvRatio", 3.0487)

    dataset.jobQueue = "gsa"

    dataset
  }

  def newK1gDataset(projectName: String, chr20: Boolean) = {
    val project = PipelineTest.createHg19Project(projectName, chr20)
    var samples = List.empty[PipelineSample]
    for (k1gBam <- PipelineTest.k1gBams)
      samples :+= PipelineTest.createK1gSample(projectName, k1gBam)
    new PipelineDataset(PipelineTest.createPipeline(project, samples))
  }

  @DataProvider(name="datasets")//, parallel=true)
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets")
  def testFullCallingPipeline(dataset: PipelineDataset) {
    val projectName = dataset.pipeline.getProject.getName
    val testName = "FullCallingPipeline-" + projectName
    val yamlFile = writeYaml(testName, dataset.pipeline)

    // Run the pipeline with the expected inputs.
    val pipelineCommand =
      "-retry 1 -S scala/qscript/playground/FullCallingPipeline.q -jobProject %s -Y %s"
        .format(projectName, yamlFile)

    val pipelineSpec = new PipelineTestSpec
    pipelineSpec.name = testName
    pipelineSpec.args = pipelineCommand
    pipelineSpec.jobQueue = dataset.jobQueue

    pipelineSpec.evalSpec = new PipelineTestEvalSpec
    pipelineSpec.evalSpec.evalReport = projectName + ".cleaned.snps_and_indels.filtered.annotated.eval"
    pipelineSpec.evalSpec.validations = dataset.validations

    PipelineTest.executeTest(pipelineSpec)
  }

  private def writeYaml(testName: String, pipeline: Pipeline) = {
    val runDir = PipelineTest.runDir(testName)
    val yamlFile = new File(runDir, pipeline.getProject.getName + ".yaml")
    yamlFile.getParentFile.mkdirs
    YamlUtils.dump(pipeline, yamlFile)
    yamlFile
  }

  class PipelineDataset(var pipeline: Pipeline = null,
                        var validations: List[PipelineValidation[_]] = Nil,
                        var jobQueue: String = null) {
    override def toString = pipeline.getProject.getName
  }
}

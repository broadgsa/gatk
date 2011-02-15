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
  def datasets = List(k1gChr20Dataset, k1gExomeDataset)

  val k1gChr20Dataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx_chr20", true)

    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.all.counter.nCalledLoci", 1348)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.known.counter.nCalledLoci", 1124)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.novel.counter.nCalledLoci", 224)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.all.titv.tiTvRatio", 3.6644)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.known.titv.tiTvRatio", 3.7426)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.novel.titv.tiTvRatio", 3.3077)

    dataset
  }

  val k1gExomeDataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx", false)

    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.all.counter.nCalledLoci", 50755)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.known.counter.nCalledLoci", 40894)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.novel.counter.nCalledLoci", 9861)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.all.titv.tiTvRatio", 3.2820)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.known.titv.tiTvRatio", 3.3384)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.novel.titv.tiTvRatio", 3.0630)

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

  @Test(dataProvider="datasets", enabled=false)
  def testFullCallingPipeline(dataset: PipelineDataset) = {
    val projectName = dataset.pipeline.getProject.getName
    val testName = "FullCallingPipeline-" + projectName
    val yamlFile = writeYaml(testName, dataset.pipeline)

    // Run the pipeline with the expected inputs.
    val pipelineCommand = ("-retry 1 -S scala/qscript/playground/FullCallingPipeline.q" +
            " -jobProject %s -Y %s" +
            " -tearScript %s/R/DataProcessingReport/GetTearsheetStats.R" +
            " --gatkjar %s")
            .format(projectName, yamlFile, PipelineTest.currentStingDir, PipelineTest.currentGATK)

    val pipelineSpec = new PipelineTestSpec
    pipelineSpec.name = testName
    pipelineSpec.args = pipelineCommand
    pipelineSpec.jobQueue = dataset.jobQueue

    pipelineSpec.evalSpec = new PipelineTestEvalSpec
    pipelineSpec.evalSpec.vcf = new File(PipelineTest.runDir(testName) + "SnpCalls/%s.cleaned.annotated.handfiltered.vcf".format(projectName))
    pipelineSpec.evalSpec.reference = dataset.pipeline.getProject.getReferenceFile
    pipelineSpec.evalSpec.intervals = dataset.pipeline.getProject.getIntervalList
    pipelineSpec.evalSpec.dbsnp = dataset.pipeline.getProject.getEvalDbsnp
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
                        var validations: List[PipelineValidation] = Nil,
                        var jobQueue: String = null) {
    override def toString = pipeline.getProject.getName
  }
}

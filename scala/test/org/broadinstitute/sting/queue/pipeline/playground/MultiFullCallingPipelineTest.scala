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

import collection.JavaConversions._
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.testng.annotations.{Test, DataProvider}
import org.broadinstitute.sting.queue.pipeline.{PipelineTestSpec, PipelineTest}
import java.io.{PrintWriter, File}
import org.apache.commons.io.IOUtils

class MultiFullCallingPipelineTest {
  def datasets = List(k1gChr20Dataset)

  val k1gChr20Dataset = newK1gDataset("Barcoded_1000G_WEx_chr20", true, "hour")
  val k1gExomeDataset = newK1gDataset("Barcoded_1000G_WEx", false, "gsa")

  def newK1gDataset(datasetName: String, chr20: Boolean, pipelineJobQueue: String) = {
    var dataset = new MultiPipelineDataset
    dataset.name = datasetName
    dataset.pipelineJobQueue = pipelineJobQueue
    for (k1gBam <- PipelineTest.k1gBams) {
      val project = PipelineTest.createHg19Project("SingleSample_" + k1gBam.sampleId, chr20)
      val sample = PipelineTest.createK1gSample("Sample", k1gBam)
      dataset.samplePipelines :+= PipelineTest.createPipeline(project, List(sample))
    }
    dataset
  }

  @DataProvider(name="datasets")//, parallel=true)
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets", enabled=false)
  def testMultiFullCallingPipeline(dataset: MultiPipelineDataset) = {
    val projectName = dataset.name
    val testName = "MultiFullCallingPipeline-" + projectName

    var yamlFiles = List.empty[File]
    for (samplePipeline <- dataset.samplePipelines)
      yamlFiles :+= writeYaml(testName, samplePipeline)

    val yamlList = writeYamlList(testName, yamlFiles)

    // Run the pipeline with the expected inputs.
    val pipelineCommand = ("-retry 1 -BS 3 -PP 100 -S scala/qscript/playground/MultiFullCallingPipeline.scala" +
            " -jobProject %s -YL %s -PJQ %s -stingHome %s")
            .format(projectName, yamlList, dataset.pipelineJobQueue, PipelineTest.currentStingDir)

    val pipelineSpec = new PipelineTestSpec
    pipelineSpec.name = testName
    pipelineSpec.args = pipelineCommand
    pipelineSpec.jobQueue = "gsa"

    PipelineTest.executeTest(pipelineSpec)
  }

  private def writeYaml(testName: String, pipeline: Pipeline) = {
    val runDir = PipelineTest.runDir(testName)
    val yamlFile = new File(runDir, pipeline.getProject.getName + "/" + pipeline.getProject.getName + ".yaml").getAbsoluteFile
    yamlFile.getParentFile.mkdirs
    YamlUtils.dump(pipeline, yamlFile)
    yamlFile
  }

  private def writeYamlList(testName: String, yamlFiles: List[File]) = {
    val runDir = PipelineTest.runDir(testName)
    val yamlList = new File(runDir, testName + "_yamls.list").getAbsoluteFile
    yamlList.getParentFile.mkdirs
    val writer = new PrintWriter(yamlList)
    try {
      for (yamlFile <- yamlFiles)
        writer.println(yamlFile.toString)
    } finally {
      IOUtils.closeQuietly(writer)
    }
    yamlList
  }

  class MultiPipelineDataset (var name: String = null,
                              var pipelineJobQueue: String = null,
                              var samplePipelines: List[Pipeline] = Nil) {
    override def toString = name
  }
}

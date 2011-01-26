package org.broadinstitute.sting.queue.pipeline.playground

import org.testng.annotations.{DataProvider, Test}
import collection.JavaConversions._
import java.io.File
import org.broadinstitute.sting.datasources.pipeline.{PipelineSample, PipelineProject, Pipeline}
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.BaseTest
import org.broadinstitute.sting.queue.pipeline._

class FullCallingPipelineTest {
  def datasets = List(k1gChr20Dataset, k1gExomeDataset)

  // In fullCallingPipeline.q VariantEval is always compared against 129.
  // Until the newvarianteval is finalized which will allow java import of the prior results,
  // we re-run VariantEval to validate the run, and replicate that behavior here.
  private final val variantEvalDbsnpFile = new File(BaseTest.b37dbSNP129)

  val k1gChr20Dataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx_chr20")
    dataset.pipeline.getProject.setIntervalList(new File(BaseTest.GATKDataLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.chr20.interval_list"))

    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.all.counter.nCalledLoci", 1348)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.known.counter.nCalledLoci", 1124)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.novel.counter.nCalledLoci", 224)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.all.titv.tiTvRatio", 3.6644)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.known.titv.tiTvRatio", 3.7426)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.novel.titv.tiTvRatio", 3.3077)

    dataset.jobQueue = "hour"

    dataset
  }

  val k1gExomeDataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx")
    dataset.pipeline.getProject.setIntervalList(new File(BaseTest.GATKDataLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"))

    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.all.counter.nCalledLoci", 50755)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.known.counter.nCalledLoci", 40894)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.novel.counter.nCalledLoci", 9861)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.all.titv.tiTvRatio", 3.2820)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.known.titv.tiTvRatio", 3.3384)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.novel.titv.tiTvRatio", 3.0630)

    dataset.jobQueue = "gsa"

    dataset
  }

  def newK1gDataset(projectName: String) = {
    val project = new PipelineProject
    project.setName(projectName)
    project.setReferenceFile(new File(BaseTest.hg19Reference))
    project.setDbsnpFile(new File(BaseTest.b37dbSNP132))
    project.setRefseqTable(new File(BaseTest.hg19Refseq))

    val squid = "C426"
    val ids = List(
      "NA19651","NA19655","NA19669","NA19834","HG01440",
      "NA12342","NA12748","NA19649","NA19652","NA19654")
    var samples = List.empty[PipelineSample]
    for (id <- ids) {
      val sample = new PipelineSample
      sample.setId(projectName + "_" + id)
      sample.setBamFiles(Map("recalibrated" -> new File("/seq/picard_aggregation/%1$s/%2$s/v6/%2$s.bam".format(squid,id))))
      sample.setTags(Map("SQUIDProject" -> squid, "CollaboratorID" -> id))
      samples :+= sample
    }

    val pipeline = new Pipeline
    pipeline.setProject(project)
    pipeline.setSamples(samples)

    val dataset = new PipelineDataset
    dataset.pipeline = pipeline

    dataset
  }

  @DataProvider(name="datasets")//, parallel=true)
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets", enabled=false)
  def testFullCallingPipeline(dataset: PipelineDataset) = {
    val projectName = dataset.pipeline.getProject.getName
    val testName = "fullCallingPipeline-" + projectName
    val yamlFile = writeTempYaml(dataset.pipeline)
    var cleanType = "cleaned"

    // Run the pipeline with the expected inputs.
    var pipelineCommand = ("-retry 1 -S scala/qscript/playground/fullCallingPipeline.q" +
            " -jobProject %s -Y %s" +
            " -tearScript %s/R/DataProcessingReport/GetTearsheetStats.R" +
            " --gatkjar %s/dist/GenomeAnalysisTK.jar")
            .format(projectName, yamlFile, PipelineTest.currentDir, PipelineTest.currentDir)

    if (!dataset.runIndelRealigner) {
      pipelineCommand += " -skipCleaning"
      cleanType = "uncleaned"
    }
    
    if (dataset.jobQueue != null)
      pipelineCommand += " -jobQueue " + dataset.jobQueue

    val pipelineSpec = new PipelineTestSpec
    pipelineSpec.args = pipelineCommand

    pipelineSpec.evalSpec = new PipelineTestEvalSpec
    pipelineSpec.evalSpec.vcf = new File(PipelineTest.runDir(testName) + "SnpCalls/%s.%s.annotated.handfiltered.vcf".format(projectName, cleanType))
    pipelineSpec.evalSpec.reference = dataset.pipeline.getProject.getReferenceFile
    pipelineSpec.evalSpec.intervals = dataset.pipeline.getProject.getIntervalList
    pipelineSpec.evalSpec.dbsnp = variantEvalDbsnpFile
    pipelineSpec.evalSpec.validations = dataset.validations

    // Run the test, at least checking if the command compiles
    PipelineTest.executeTest(testName, pipelineSpec)
  }

  class PipelineDataset(
          var pipeline: Pipeline = null,
          var validations: List[PipelineValidation] = Nil,
          var jobQueue: String = null,
          var runIndelRealigner: Boolean = false) {
    override def toString = pipeline.getProject.getName
  }

  private def writeTempYaml(pipeline: Pipeline) = {
    val tempFile = BaseTest.createTempFile(pipeline.getProject.getName + "-", ".yaml")
    YamlUtils.dump(pipeline, tempFile)
    tempFile
  }
}

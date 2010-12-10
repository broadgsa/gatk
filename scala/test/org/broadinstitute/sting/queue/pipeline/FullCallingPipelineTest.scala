package org.broadinstitute.sting.queue.pipeline

import org.testng.annotations.{DataProvider, Test}
import collection.JavaConversions._
import java.io.File
import org.broadinstitute.sting.datasources.pipeline.{PipelineSample, PipelineProject, Pipeline}
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.{WalkerTest, BaseTest}
import java.text.SimpleDateFormat
import java.util.Date

class FullCallingPipelineTest extends BaseTest {
  def datasets = List(k1gChr20Dataset, k1gExomeDataset)

  private final val validationReportsDataLocation = "/humgen/gsa-hpprojects/GATK/validationreports/submitted/"

  val k1gChr20Dataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx_chr20")
    dataset.pipeline.getProject.setIntervalList(new File(BaseTest.GATKDataLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.chr20.interval_list"))

    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.all.counter.nCalledLoci", 1359)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.known.counter.nCalledLoci", 1134)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.novel.counter.nCalledLoci", 225)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.all.titv.tiTvRatio", 3.6701)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.known.titv.tiTvRatio", 3.7647)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.novel.titv.tiTvRatio", 3.2453)

    dataset.jobQueue = "hour"

    dataset
  }

  val k1gExomeDataset = {
    val dataset = newK1gDataset("Barcoded_1000G_WEx")
    dataset.pipeline.getProject.setIntervalList(new File(BaseTest.GATKDataLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"))

    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.all.counter.nCalledLoci", 51130)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.known.counter.nCalledLoci", 41042)
    dataset.validations :+= new IntegerValidation("eval.dbsnp.all.called.novel.counter.nCalledLoci", 10088)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.all.titv.tiTvRatio", 3.2598)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.known.titv.tiTvRatio", 3.3307)
    dataset.validations :+= new DoubleValidation("eval.dbsnp.all.called.novel.titv.tiTvRatio", 2.9937)

    dataset.jobQueue = "gsa"
    dataset.bigMemQueue = "gsa"

    dataset
  }

  def newK1gDataset(projectName: String) = {
    val project = new PipelineProject
    project.setName(projectName)
    project.setReferenceFile(new File(BaseTest.hg19Reference))
    project.setDbsnpFile(new File(BaseTest.b37dbSNP129))

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
    dataset.refseq = BaseTest.hg19Refseq
    dataset.targetTiTv = "3.0"

    dataset
  }

  @DataProvider(name="datasets")//, parallel=true)
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets")
  def testFullCallingPipeline(dataset: PipelineDataset) = {
    val projectName = dataset.pipeline.getProject.getName
    val testName = "fullCallingPipeline-" + projectName
    val yamlFile = writeTempYaml(dataset.pipeline)
    var cleanType = "cleaned"

    // Run the pipeline with the expected inputs.
    var pipelineCommand = ("-retry 1 -S scala/qscript/fullCallingPipeline.q" +
            " -jobProject %s -Y %s -refseqTable %s -titv %s" +
            " --gatkjar %s/dist/GenomeAnalysisTK.jar")
            .format(projectName, yamlFile, dataset.refseq, dataset.targetTiTv, new File(".").getAbsolutePath)

    if (!dataset.runIndelRealigner) {
      pipelineCommand += " -skipCleaning"
      cleanType = "uncleaned"
    }
    
    if (dataset.jobQueue != null)
      pipelineCommand += " -jobQueue " + dataset.jobQueue

    if (dataset.bigMemQueue != null)
      pipelineCommand += " -bigMemQueue " + dataset.bigMemQueue

    // Run the test, at least checking if the command compiles
    PipelineTest.executeTest(testName, pipelineCommand, null)

    // If actually running, evaluate the output validating the expressions.
    if (PipelineTest.run) {
      // path where the pipeline should have output the handfiltered vcf
      val handFilteredVcf = PipelineTest.runDir(testName) + "SnpCalls/%s.%s.annotated.handfiltered.vcf".format(projectName, cleanType)

      // eval modules to record in the validation directory
      val evalModules = List("CompOverlap", "CountFunctionalClasses", "CountVariants", "SimpleMetricsBySample", "TiTvVariantEvaluator")

      // write the report to the shared validation data location
      val formatter = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss")
      val reportLocation = "%s%s/validation.%s.eval".format(validationReportsDataLocation, testName, formatter.format(new Date))
      new File(reportLocation).getParentFile.mkdirs

      // Run variant eval generating the report and validating the pipeline vcfs.
      var walkerCommand = ("-T VariantEval -R %s -D %s -B:eval,VCF %s" +
              " -E %s -reportType R -reportLocation %s -L %s")
              .format(
        dataset.pipeline.getProject.getReferenceFile, dataset.pipeline.getProject.getDbsnpFile, handFilteredVcf,
        evalModules.mkString(" -E "), reportLocation, dataset.pipeline.getProject.getIntervalList)

      for (validation <- dataset.validations) {
        walkerCommand += " -summary %s".format(validation.metric)
        walkerCommand += " -validate '%1$s >= %2$s' -validate '%1$s <= %3$s'".format(
          validation.metric, validation.min, validation.max)
      }

      WalkerTest.executeTest("fullCallingPipelineValidate-" + projectName, walkerCommand, null)
    }
  }

  class PipelineDataset(
          var pipeline: Pipeline = null,
          var refseq: String = null,
          var targetTiTv: String = null,
          var validations: List[PipelineValidation] = Nil,
          var jobQueue: String = null,
          var bigMemQueue: String = null,
          var runIndelRealigner: Boolean = false) {
    override def toString = pipeline.getProject.getName
  }

  class PipelineValidation(val metric: String, val min: String, val max: String) {
  }

  class IntegerValidation(metric: String, target: Int)
          extends PipelineValidation(metric,
            (target * .99).floor.toInt.toString, (target * 1.01).ceil.toInt.toString) {
  }

  class DoubleValidation(metric: String, target: Double)
          extends PipelineValidation(metric,
            "%.2f".format((target * 99).floor / 100), "%.2f".format((target * 101).ceil / 100)) {
  }

  private def writeTempYaml(pipeline: Pipeline) = {
    val tempFile = File.createTempFile(pipeline.getProject.getName + "-", ".yaml")
    tempFile.deleteOnExit
    YamlUtils.dump(pipeline, tempFile)
    tempFile
  }
}

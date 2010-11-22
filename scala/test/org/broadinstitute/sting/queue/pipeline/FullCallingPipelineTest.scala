package org.broadinstitute.sting.queue.pipeline

import org.testng.annotations.{DataProvider, Test}
import collection.JavaConversions._
import java.io.File
import org.broadinstitute.sting.datasources.pipeline.{PipelineSample, PipelineProject, Pipeline}
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.queue.PipelineTest
import org.broadinstitute.sting.{WalkerTest, BaseTest}
import java.text.SimpleDateFormat
import java.util.Date

class FullCallingPipelineTest extends BaseTest {
  def datasets = List(k1gChr20Dataset)

  private final val validationReportsDataLocation = "/humgen/gsa-hpprojects/GATK/validationreports/submitted/"

  val k1gChr20Dataset = {
    val project = new PipelineProject
    project.setName("Barcoded_1000G_WEx_chr20")
    project.setReferenceFile(new File(BaseTest.hg19Reference))
    project.setDbsnpFile(new File(BaseTest.b37dbSNP129))
    project.setIntervalList(new File(BaseTest.GATKDataLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.chr20.interval_list"))

    val squid = "C426"
    val ids = List(
      "NA19651","NA19655","NA19669","NA19834","HG01440",
      "NA12342","NA12748","NA19649","NA19652","NA19654")
    var samples = List.empty[PipelineSample]
    for (id <- ids) {
      val sample = new PipelineSample
      sample.setId(project.getName + "_" + id)
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
    dataset.jobQueue = "hour"

    dataset.validations :+= new PipelineValidation("evalHandFiltered.dbsnp.all.called.all.counter.nCalledLoci", "1390", "1420")
    dataset.validations :+= new PipelineValidation("evalHandFiltered.dbsnp.all.called.all.titv.tiTvRatio", "3.52", "3.60")
    dataset.validations :+= new PipelineValidation("evalHandFiltered.dbsnp.all.called.known.titv.tiTvRatio", "3.71", "3.80")
    dataset.validations :+= new PipelineValidation("evalHandFiltered.dbsnp.all.called.novel.titv.tiTvRatio", "2.79", "2.86")

    dataset
  }

  @DataProvider(name="datasets")
  final def convertDatasets: Array[Array[AnyRef]] =
    datasets.map(dataset => Array(dataset.asInstanceOf[AnyRef])).toArray

  @Test(dataProvider="datasets")
  def testFullCallingPipeline(dataset: PipelineDataset) = {
    val projectName = dataset.pipeline.getProject.getName
    val testName = "fullCallingPipeline-" + projectName
    val yamlFile = writeTempYaml(dataset.pipeline)

    // Run the pipeline with the expected inputs.
    var pipelineCommand = ("-jobProject %s -S scala/qscript/fullCallingPipeline.q -Y %s" +
            " -refseqTable %s" +
            " --gatkjar %s/dist/GenomeAnalysisTK.jar -titv %s -skipCleaning")
            .format(projectName, yamlFile, dataset.refseq, new File(".").getCanonicalPath, dataset.targetTiTv)

    if (dataset.jobQueue != null)
      pipelineCommand += " -jobQueue " + dataset.jobQueue

    // Run the test, at least checking if the command compiles
    PipelineTest.executeTest(testName, pipelineCommand, null)

    // If actually running, evaluate the output validating the expressions.
    if (PipelineTest.run) {
      // path where the pipeline should have output the uncleaned handfiltered vcf
      val handFilteredVcf = PipelineTest.runDir(testName) + "SnpCalls/%s.uncleaned.annotated.handfiltered.vcf".format(projectName)

      // path where the pipeline should have outout the indel masked vcf
      val optimizedVcf = PipelineTest.runDir(testName) + "SnpCalls/%s.uncleaned.annotated.indel.masked.recalibrated.tranched.vcf".format(projectName)

      // eval modules to record in the validation directory
      val evalModules = List("CountFunctionalClasses", "CompOverlap", "CountVariants", "TiTvVariantEvaluator")

      // write the report to the shared validation data location
      val formatter = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss")
      val reportLocation = "%s/%s/validation.%s.eval".format(validationReportsDataLocation, testName, formatter.format(new Date))
      new File(reportLocation).getParentFile.mkdirs

      // Run variant eval generating the report and validating the pipeline vcfs.
      var walkerCommand = ("-T VariantEval -R %s -D %s -B:evalOptimized,VCF %s -B:evalHandFiltered,VCF %s" +
              " -E %s -reportType R -reportLocation %s -L %s")
              .format(
        dataset.pipeline.getProject.getReferenceFile, dataset.pipeline.getProject.getDbsnpFile, optimizedVcf, handFilteredVcf,
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
          var jobQueue: String = null) {
    override def toString = pipeline.getProject.getName
  }

  class PipelineValidation(
          var metric: String = null,
          var min: String = null,
          var max: String = null) {
  }
  
  private def writeTempYaml(pipeline: Pipeline) = {
    val tempFile = File.createTempFile(pipeline.getProject.getName + "-", ".yaml")
    tempFile.deleteOnExit
    YamlUtils.dump(pipeline, tempFile)
    tempFile
  }
}

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

package org.broadinstitute.sting.queue.pipeline

import collection.JavaConversions._
import org.broadinstitute.sting.utils.Utils
import org.testng.Assert
import org.broadinstitute.sting.commandline.CommandLineProgram
import java.util.Date
import java.text.SimpleDateFormat
import org.broadinstitute.sting.BaseTest
import org.broadinstitute.sting.queue.QCommandLine
import org.broadinstitute.sting.datasources.pipeline.{Pipeline, PipelineProject, PipelineSample}
import org.broadinstitute.sting.utils.broad.PicardAggregationUtils
import org.broadinstitute.sting.queue.util.{Logging, ProcessController}
import java.io.{FileNotFoundException, File}
import org.broadinstitute.sting.gatk.report.GATKReportParser
import org.apache.commons.io.FileUtils
import org.broadinstitute.sting.queue.engine.CommandLinePluginManager

object PipelineTest extends BaseTest with Logging {

  case class K1gBam(squidId: String, sampleId: String, version: Int)

  /** 1000G BAMs used for validation */
  val k1gBams = List(
    new K1gBam("C474", "NA19651", 2),
    new K1gBam("C474", "NA19655", 2),
    new K1gBam("C474", "NA19669", 2),
    new K1gBam("C454", "NA19834", 2),
    new K1gBam("C460", "HG01440", 2),
    new K1gBam("C456", "NA12342", 2),
    new K1gBam("C456", "NA12748", 2),
    new K1gBam("C474", "NA19649", 2),
    new K1gBam("C474", "NA19652", 2),
    new K1gBam("C474", "NA19654", 2))

  validateK1gBams()

  private val validationReportsDataLocation = "/humgen/gsa-hpprojects/GATK/validationreports/submitted/"

  val run = System.getProperty("pipeline.run") == "run"

  private val jobRunners = {
    val commandLinePluginManager = new CommandLinePluginManager
    commandLinePluginManager.getPlugins.map(commandLinePluginManager.getName(_)).filterNot(_ == "Shell")
  }

  /**
   * Returns the top level output path to this test.
   * @param testName The name of the test passed to PipelineTest.executeTest()
   * @param jobRunner The name of the job manager to run the jobs.
   * @return the top level output path to this test.
   */
  def testDir(testName: String, jobRunner: String) = "pipelinetests/%s/%s/".format(testName, jobRunner)

  /**
   * Returns the directory where relative output files will be written for this test.
   * @param testName The name of the test passed to PipelineTest.executeTest()
   * @param jobRunner The name of the job manager to run the jobs.
   * @return the directory where relative output files will be written for this test.
   */
  private def runDir(testName: String, jobRunner: String) = testDir(testName, jobRunner) + "run/"

  /**
   * Returns the directory where temp files will be written for this test.
   * @param testName The name of the test passed to PipelineTest.executeTest()
   * @param jobRunner The name of the job manager to run the jobs.
   * @return the directory where temp files will be written for this test.
   */
  private def tempDir(testName: String, jobRunner: String) = testDir(testName, jobRunner) + "temp/"

  /**
   * Creates a new pipeline from a project.
   * @param project Pipeline project info.
   * @param samples List of samples.
   * @return a new pipeline project.
   */
  def createPipeline(project: PipelineProject, samples: List[PipelineSample]) = {
    val pipeline = new Pipeline
    pipeline.setProject(project)
    pipeline.setSamples(samples)
    pipeline
  }

  /**
   * Creates a new pipeline project for hg19 with b37 132 dbsnp for genotyping, and b37 129 dbsnp for eval.
   * @param projectName Name of the project.
   * @param intervals The intervals file to use.
   * @return a new pipeline project.
   */
  def createHg19Project(projectName: String, intervals: String) = {
    val project = new PipelineProject
    project.setName(projectName)
    project.setReferenceFile(new File(BaseTest.hg19Reference))
    project.setGenotypeDbsnp(new File(BaseTest.b37dbSNP132))
    project.setEvalDbsnp(new File(BaseTest.b37dbSNP129))
    project.setRefseqTable(new File(BaseTest.hg19Refseq))
    project.setIntervalList(new File(intervals))
    project
  }

  /**
   * Creates a 1000G pipeline sample from one of the bams.
   * @param idPrefix Text to prepend to the sample name.
   * @param k1gBam bam to create the sample for.
   * @return the created pipeline sample.
   */
  def createK1gSample(idPrefix: String, k1gBam: K1gBam) = {
    val sample = new PipelineSample
    sample.setId(idPrefix + "_" + k1gBam.sampleId)
    sample.setBamFiles(Map("cleaned" -> getPicardBam(k1gBam)))
    sample
  }

  /**
   * Runs the pipelineTest.
   * @param pipelineTest test to run.
   */
  def executeTest(pipelineTest: PipelineTestSpec) {
    jobRunners.foreach(executeTest(pipelineTest, _))
  }
  
  /**
   * Runs the pipelineTest.
   * @param pipelineTest test to run.
   * @param jobRunner The name of the job manager to run the jobs.
   */
  def executeTest(pipelineTest: PipelineTestSpec, jobRunner: String) {
    val name = pipelineTest.name
    if (name == null)
      Assert.fail("PipelineTestSpec.name is null")
    println(Utils.dupString('-', 80));
    executeTest(name, pipelineTest.args, pipelineTest.jobQueue, pipelineTest.expectedException, jobRunner)
    if (run) {
      assertMatchingMD5s(name, pipelineTest.fileMD5s.map{case (file, md5) => new File(runDir(name, jobRunner), file) -> md5}, pipelineTest.parameterize)
      if (pipelineTest.evalSpec != null)
        validateEval(name, pipelineTest.evalSpec, jobRunner)
      println("  => %s PASSED (%s)".format(name, jobRunner))
    }
    else
      println("  => %s PASSED DRY RUN (%s)".format(name, jobRunner))
  }

  private def assertMatchingMD5s(name: String, fileMD5s: Traversable[(File, String)], parameterize: Boolean) {
    var failed = 0
    for ((file, expectedMD5) <- fileMD5s) {
      val calculatedMD5 = BaseTest.testFileMD5(name, file, expectedMD5, parameterize)
      if (!parameterize && expectedMD5 != "" && expectedMD5 != calculatedMD5)
        failed += 1
    }
    if (failed > 0)
      Assert.fail("%d of %d MD5s did not match".format(failed, fileMD5s.size))
  }

  private def validateEval(name: String, evalSpec: PipelineTestEvalSpec, jobRunner: String) {
    // write the report to the shared validation data location
    val formatter = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss")
    val reportLocation = "%s%s/%s/validation.%s.eval".format(validationReportsDataLocation, jobRunner, name, formatter.format(new Date))
    val report = new File(reportLocation)

    FileUtils.copyFile(new File(runDir(name, jobRunner) + evalSpec.evalReport), report);

    val parser = new GATKReportParser
    parser.parse(report)

    var allInRange = true

    println()
    println(name + " validation values:")
    println("    value (min,target,max) table key metric")
    for (validation <- evalSpec.validations) {
      val value = parser.getValue(validation.table, validation.key, validation.metric)
      val inRange = if (value == null) false else validation.inRange(value)
      val flag = if (!inRange) "*" else " "
      println("  %s %s (%s,%s,%s) %s %s %s".format(flag, value, validation.min, validation.target, validation.max, validation.table, validation.key, validation.metric))
      allInRange &= inRange
    }

    if (!allInRange)
      Assert.fail("Eval outside of expected range")
  }

  /**
   * execute the test
   * @param name the name of the test
   * @param args the argument list
   * @param jobQueue the queue to run the job on.  Defaults to hour if jobQueue is null.
   * @param expectedException the expected exception or null if no exception is expected.
   * @param jobRunner The name of the job manager to run the jobs.
   */
  private def executeTest(name: String, args: String, jobQueue: String, expectedException: Class[_], jobRunner: String) {
    var command = Utils.escapeExpressions(args)

    // add the logging level to each of the integration test commands

    command = Utils.appendArray(command, "-jobRunner", jobRunner,
      "-tempDir", tempDir(name, jobRunner), "-runDir", runDir(name, jobRunner))

    if (jobQueue != null)
      command = Utils.appendArray(command, "-jobQueue", jobQueue)

    if (run)
      command = Utils.appendArray(command, "-run")

    // run the executable
    var gotAnException = false

    val instance = new QCommandLine
    runningCommandLines += instance
    try {
      println("Executing test %s with Queue arguments: %s".format(name, Utils.join(" ",command)))
      CommandLineProgram.start(instance, command)
    } catch {
      case e =>
        gotAnException = true
        if (expectedException != null) {
          // we expect an exception
          println("Wanted exception %s, saw %s".format(expectedException, e.getClass))
          if (expectedException.isInstance(e)) {
            // it's the type we expected
            println(String.format("  => %s PASSED (%s)", name, jobRunner))
          } else {
            e.printStackTrace()
            Assert.fail("Test %s expected exception %s but got %s instead (%s)".format(
              name, expectedException, e.getClass, jobRunner))
          }
        } else {
          // we didn't expect an exception but we got one :-(
          throw new RuntimeException(e)
        }
    } finally {
      instance.shutdown()
      runningCommandLines -= instance
    }

    // catch failures from the integration test
    if (expectedException != null) {
      if (!gotAnException)
      // we expected an exception but didn't see it
        Assert.fail("Test %s expected exception %s but none was thrown (%s)".format(name, expectedException.toString, jobRunner))
    } else {
      if (CommandLineProgram.result != 0)
        throw new RuntimeException("Error running Queue with arguments: " + args)
    }
  }

  /**
   * Throws an exception if any of the 1000G bams do not exist and warns if they are out of date.
   */
  private def validateK1gBams() {
    var missingBams = List.empty[File]
    for (k1gBam <- k1gBams) {
      val latest = getLatestVersion(k1gBam)
      val bam = getPicardBam(k1gBam)
      if (k1gBam.version != latest)
        logger.warn("1000G bam is not the latest version %d: %s".format(latest, k1gBam))
      if (!bam.exists)
        missingBams :+= bam
    }
    if (missingBams.size > 0) {
      val nl = "%n".format()
      throw new FileNotFoundException("The following 1000G bam files are missing.%n%s".format(missingBams.mkString(nl)))
    }
  }

  private def getPicardBam(k1gBam: K1gBam): File =
    new File(PicardAggregationUtils.getSampleBam(k1gBam.squidId, k1gBam.sampleId, k1gBam.version))

  private def getLatestVersion(k1gBam: K1gBam): Int =
    PicardAggregationUtils.getLatestVersion(k1gBam.squidId, k1gBam.sampleId, k1gBam.version)

  private var runningCommandLines = Set.empty[QCommandLine]

  Runtime.getRuntime.addShutdownHook(new Thread {
    /** Cleanup as the JVM shuts down. */
    override def run() {
      try {
        ProcessController.shutdown()
      } catch {
        case _ => /*ignore */
      }
      runningCommandLines.foreach(commandLine =>
        try {
          commandLine.shutdown()
        } catch {
          case _ => /* ignore */
        })
    }
  })
}

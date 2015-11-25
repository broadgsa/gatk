/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.pipeline

import org.broadinstitute.gatk.utils.Utils
import org.testng.Assert
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram
import java.util.Date
import java.text.SimpleDateFormat
import org.broadinstitute.gatk.utils.BaseTest
import org.broadinstitute.gatk.utils.MD5DB
import org.broadinstitute.gatk.queue.{QScript, QCommandLine}
import org.broadinstitute.gatk.queue.util.Logging
import java.io.{FilenameFilter, File}
import org.apache.commons.io.FileUtils
import org.apache.commons.io.filefilter.WildcardFileFilter
import org.broadinstitute.gatk.utils.report.GATKReport

object QueueTest extends BaseTest with Logging {

  private final val qscriptsSrcDir = "src/main/qscripts/"
  private final val qscriptsPackageDir = qscriptsSrcDir + "org/broadinstitute/gatk/queue/qscripts/"
  final val publicQScriptsPackageDir = BaseTest.gatkDirectory + "public/gatk-queue-extensions-public/" + qscriptsPackageDir
  final val protectedQScriptsPackageDir = BaseTest.gatkDirectory + "protected/gatk-queue-extensions-distribution/" + qscriptsPackageDir
  final val privateQScriptsPackageDir = BaseTest.gatkDirectory + "private/gatk-queue-extensions-internal/" + qscriptsPackageDir

  private val validationReportsDataLocation = "/humgen/gsa-hpprojects/GATK/validationreports/submitted/"
  private val md5DB = new MD5DB()

  /**
   * All the job runners configured to run QueueTests at The Broad.
   */
  final val allJobRunners = Seq("GridEngine", "Shell", "ParallelShell")

  /**
   * The default job runners to run.
   */
  final val defaultJobRunners = Seq("GridEngine")

  /**
   * Returns the top level output path to this test.
   * @param testName The name of the test passed to QueueTest.executeTest()
   * @param jobRunner The name of the job manager to run the jobs.
   * @return the top level output path to this test.
   */
  def testDir(testName: String, jobRunner: String) = "queuetests/%s/%s/".format(testName, jobRunner)

  /**
   * Returns the directory where relative output files will be written for this test.
   * @param testName The name of the test passed to QueueTest.executeTest()
   * @param jobRunner The name of the job manager to run the jobs.
   * @return the directory where relative output files will be written for this test.
   */
  private def runDir(testName: String, jobRunner: String) = testDir(testName, jobRunner) + "run/"

  /**
   * Returns the directory where temp files will be written for this test.
   * @param testName The name of the test passed to QueueTest.executeTest()
   * @param jobRunner The name of the job manager to run the jobs.
   * @return the directory where temp files will be written for this test.
   */
  private def tempDir(testName: String, jobRunner: String) = testDir(testName, jobRunner) + "temp/"

  /**
   * Runs the queueTest.
   * @param queueTest test to run.
   */
  def executeTest(queueTest: QueueTestSpec) {
    var jobRunners = queueTest.jobRunners
    if (jobRunners == null)
      jobRunners = defaultJobRunners
    jobRunners.foreach(executeTest(queueTest, _))
  }

  /**
   * Runs the queueTest.
   * @param queueTest test to run.
   * @param jobRunner The name of the job manager to run the jobs.
   */
  def executeTest(queueTest: QueueTestSpec, jobRunner: String) {
    // Reset the order of functions added to the graph.
    QScript.resetAddOrder()

    val name = queueTest.name
    if (name == null)
      Assert.fail("QueueTestSpec.name is null")
    println(Utils.dupString('-', 80))
    executeTest(name, queueTest.args, queueTest.jobQueue, queueTest.expectedException, jobRunner)
    if (BaseTest.queueTestRunModeIsSet) {
      assertMatchingMD5s(name, queueTest.fileMD5s.map{case (file, md5) => new File(runDir(name, jobRunner), file) -> md5}, queueTest.parameterize)
      if (queueTest.evalSpec != null)
        validateEval(name, queueTest.evalSpec, jobRunner)
      for (path <- queueTest.expectedFilePaths)
        assertPathExists(runDir(name, jobRunner), path)
      for (path <- queueTest.unexpectedFilePaths)
        assertPathDoesNotExist(runDir(name, jobRunner), path)
      println("  => %s PASSED (%s)".format(name, jobRunner))
    }
    else
      println("  => %s PASSED DRY RUN (%s)".format(name, jobRunner))
  }

  private def assertMatchingMD5s(name: String, fileMD5s: Traversable[(File, String)], parameterize: Boolean) {
    var failed = 0
    for ((file, expectedMD5) <- fileMD5s) {
      val calculatedMD5 = md5DB.testFileMD5(name, "", file, expectedMD5, parameterize).actualMD5
      if (!parameterize && expectedMD5 != "" && expectedMD5 != calculatedMD5)
        failed += 1
    }
    if (failed > 0)
      Assert.fail("%d of %d MD5s did not match".format(failed, fileMD5s.size))
  }

  private def validateEval(name: String, evalSpec: QueueTestEvalSpec, jobRunner: String) {
    // write the report to the shared validation data location
    val formatter = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss")
    val reportLocation = "%s%s/%s/validation.%s.eval".format(validationReportsDataLocation, jobRunner, name, formatter.format(new Date))
    val reportFile = new File(reportLocation)

    FileUtils.copyFile(new File(runDir(name, jobRunner) + evalSpec.evalReport), reportFile)

    val report = new GATKReport(reportFile)

    var allInRange = true

    println()
    println(name + " validation values:")
    println("    value (min,target,max) table key metric")
    for (validation <- evalSpec.validations) {
      val table = report.getTable(validation.table)
      val key = table.findRowByData(validation.table +: validation.key.split('.') : _*)
      val value = String.valueOf(table.get(key, validation.metric))
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

    if (BaseTest.queueTestRunModeIsSet)
      command = Utils.appendArray(command, "-run")

    // run the executable
    var gotAnException = false

    val instance = new QCommandLine
    runningCommandLines += instance
    try {
      println("Executing test %s with Queue arguments: %s".format(name, Utils.join(" ",command)))
      CommandLineProgram.start(instance, command)
    } catch {
      case e: Exception =>
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

  private def assertPathExists(runDir: String, path: String) {
    val orig = new File(runDir, path)
    var dir = orig.getParentFile
    if (dir == null)
      dir = new File(".")
    Assert.assertTrue(dir.exists, "Missing directory: " + dir.getAbsolutePath)
    val filter: FilenameFilter = new WildcardFileFilter(orig.getName)
    Assert.assertNotEquals(dir.listFiles(filter).length, 0, "Missing file: " + orig.getAbsolutePath)
  }

  private def assertPathDoesNotExist(runDir: String, path: String) {
    val orig = new File(runDir, path)
    var dir = orig.getParentFile
    if (dir == null)
      dir = new File(".")
    if (dir.exists) {
      val filter: FilenameFilter = new WildcardFileFilter(orig.getName)
      Assert.assertEquals(dir.listFiles(filter).length, 0,
        "Found unexpected file(s): " + dir.listFiles().map(_.getAbsolutePath).mkString(", "))
    }
  }

  private var runningCommandLines = Set.empty[QCommandLine]

  Runtime.getRuntime.addShutdownHook(new Thread {
    /** Cleanup as the JVM shuts down. */
    override def run() {
      runningCommandLines.foreach(commandLine =>
        try {
          commandLine.shutdown()
        } catch {
          case _: Throwable => /* ignore */
        })
    }
  })
}

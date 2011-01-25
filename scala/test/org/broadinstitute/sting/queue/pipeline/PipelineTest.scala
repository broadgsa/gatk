package org.broadinstitute.sting.queue.pipeline

import org.broadinstitute.sting.utils.Utils
import org.testng.Assert
import org.broadinstitute.sting.commandline.CommandLineProgram
import org.broadinstitute.sting.queue.QCommandLine
import java.io.File
import org.broadinstitute.sting.queue.util.{TextFormatUtils, ProcessController}
import java.util.Date
import java.text.SimpleDateFormat
import org.broadinstitute.sting.{WalkerTest, BaseTest}

object PipelineTest {
  private var runningCommandLines = Set.empty[QCommandLine]

  private val validationReportsDataLocation = "/humgen/gsa-hpprojects/GATK/validationreports/submitted/"

  val run = System.getProperty("pipeline.run") == "run"

  def executeTest(name: String, pipelineTest: PipelineTestSpec) {
    println(Utils.dupString('-', 80));
    executeTest(name, pipelineTest.args, pipelineTest.expectedException)
    if (run) {
      assertMatchingMD5s(name, pipelineTest.fileMD5s)
      if (pipelineTest.evalSpec != null)
        validateEval(name, pipelineTest.evalSpec)
      println("  => %s PASSED".format(name))
    }
    else
      println("  => %s PASSED DRY RUN".format(name))
  }

  private def assertMatchingMD5s(name: String, fileMD5s: Traversable[(File, String)]) {
    var failed = 0
    for ((file, expectedMD5) <- fileMD5s) {
      val calculatedMD5 = BaseTest.testFileMD5(name, file, expectedMD5, false)
      failed += 1
    }
    if (failed > 0)
      Assert.fail("%d MD5%s did not match.".format(failed, TextFormatUtils.plural(failed)))
  }

  private def validateEval(name: String, evalSpec: PipelineTestEvalSpec) {
    // write the report to the shared validation data location
    val formatter = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss")
    val reportLocation = "%s%s/validation.%s.eval".format(validationReportsDataLocation, name, formatter.format(new Date))
    new File(reportLocation).getParentFile.mkdirs

    // Run variant eval generating the report and validating the pipeline vcf.
    var walkerCommand = "-T VariantEval -R %s -B:eval,VCF %s -E %s -reportType R -reportLocation %s -L %s"
      .format(evalSpec.reference, evalSpec.vcf, evalSpec.evalModules.mkString(" -E "), reportLocation, evalSpec.intervals)

    if (evalSpec.dbsnp != null) {
      val dbsnpArg = if (evalSpec.dbsnp.getName.toLowerCase.endsWith(".vcf")) "-B:dbsnp,VCF" else "-D"
      walkerCommand += " %s %s".format(dbsnpArg, evalSpec.dbsnp)
    }

    if (evalSpec.intervals != null)
      walkerCommand += " -L %s".format(evalSpec.intervals)

    for (validation <- evalSpec.validations) {
      walkerCommand += " -summary %s".format(validation.metric)
      walkerCommand += " -validate '%1$s >= %2$s' -validate '%1$s <= %3$s'".format(
        validation.metric, validation.min, validation.max)
    }

    WalkerTest.executeTest(name + "-validate", walkerCommand, null)
  }

  /**
   * execute the test
   * @param name the name of the test
   * @param args the argument list
   * @param expectedException the expected exception or null if no exception is expected.
   */
  def executeTest(name: String, args: String, expectedException: Class[_]) {
    var command = Utils.escapeExpressions(args)

    // add the logging level to each of the integration test commands

    command = Utils.appendArray(command, "-bsub", "-l", "WARN", "-startFromScratch", "-tempDir", tempDir(name), "-runDir", runDir(name))
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
            println(String.format("  => %s PASSED", name))
          } else {
            e.printStackTrace()
            Assert.fail("Test %s expected exception %s but got %s instead".format(
              name, expectedException, e.getClass))
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
        Assert.fail("Test %s expected exception %s but none was thrown".format(name, expectedException.toString))
    } else {
      if (CommandLineProgram.result != 0)
        throw new RuntimeException("Error running the GATK with arguments: " + args)
    }
  }

  val currentDir = new File(".").getAbsolutePath
  def testDir(testName: String) = "pipelinetests/%s/".format(testName)
  def runDir(testName: String) = testDir(testName) + "run/"
  def tempDir(testName: String) = testDir(testName) + "temp/"

  Runtime.getRuntime.addShutdownHook(new Thread {
    /** Cleanup as the JVM shuts down. */
    override def run {
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

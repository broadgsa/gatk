package org.broadinstitute.sting.queue

import org.broadinstitute.sting.utils.Utils
import org.testng.Assert
import org.broadinstitute.sting.commandline.CommandLineProgram
import org.broadinstitute.sting.queue.util.ProcessController

object PipelineTest {
  private var runningCommandLines = Set.empty[QCommandLine]

  val run = System.getProperty("pipeline.run") == "run"

  /**
   * execute the test
   * @param name the name of the test
   * @param args the argument list
   * @param expectedException the expected exception or null if no exception is expected.
   */
  def executeTest(name: String, args: String, expectedException: Class[_]) = {
    var command = Utils.escapeExpressions(args)

    // add the logging level to each of the integration test commands

    command = Utils.appendArray(command, "-l", "WARN", "-startFromScratch", "-tempDir", tempDir(name), "-runDir", runDir(name))
    if (run)
      command = Utils.appendArray(command, "-run", "-bsub")

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

  def testDir(testName: String) = "pipelinetests/%s/".format(testName)
  def runDir(testName: String) = testDir(testName) + "run/"
  def tempDir(testName: String) = testDir(testName) + "temp/"

  Runtime.getRuntime.addShutdownHook(new Thread {
    /** Cleanup as the JVM shuts down. */
    override def run = {
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

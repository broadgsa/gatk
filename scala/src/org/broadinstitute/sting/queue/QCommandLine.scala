package org.broadinstitute.sting.queue

import java.io.File
import java.util.Arrays
import org.broadinstitute.sting.queue.engine.QGraph
import org.broadinstitute.sting.commandline.{ClassType, Input, Argument, CommandLineProgram}
import org.broadinstitute.sting.queue.util.{Logging, ScalaCompoundArgumentTypeDescriptor}

/**
 * Entry point of Queue.  Compiles and runs QScripts passed in to the command line.
 */
class QCommandLine extends CommandLineProgram with Logging {
  @Input(fullName="script", shortName="S", doc="QScript scala file", required=true)
  @ClassType(classOf[File])
  private var scripts = List.empty[File]

  @Argument(fullName="bsub_all_jobs", shortName="bsub", doc="Use bsub to submit jobs", required=false)
  private var bsubAllJobs = false

  @Argument(fullName="bsub_wait_jobs", shortName="bsubWait", doc="Wait for bsub submitted jobs before exiting", required=false)
  private var bsubWaitJobs = false

  @Argument(fullName="run_scripts", shortName="run", doc="Run QScripts", required=false)
  private var run = false

  @Argument(fullName="dot_graph", shortName="dot", doc="Outputs the queue graph to a .dot file.  See: http://en.wikipedia.org/wiki/DOT_language", required=false)
  private var queueDot: File = _

  /**
   * Takes the QScripts passed in, runs their script() methods, retrieves their generated
   * functions, and then builds and runs a QGraph based on the dependencies.
   */
  def execute = {
    val qGraph = new QGraph
    qGraph.dryRun = !run
    qGraph.bsubAllJobs = bsubAllJobs
    qGraph.bsubWaitJobs = bsubWaitJobs

    val scripts = qScriptManager.createScripts()
    for (script <- scripts) {
      logger.info("Scripting " + qScriptManager.getName(script.getClass.asSubclass(classOf[QScript])))
      loadArgumentsIntoObject(script)
      script.script
      script.functions.foreach(qGraph.add(_))
      logger.info("Added " + script.functions.size + " functions")
    }

    logger.info("Binding functions")
    qGraph.fillIn
    if (queueDot != null) {
      logger.info("Generating " + queueDot)
      qGraph.renderToDot(queueDot)
    }

    logger.info("Running generated graph")
    qGraph.run
    logger.info("Done")
    0
  }

  /**
   * Returns true as QScripts are located and compiled.
   * @return true
   */
  override def canAddArgumentsDynamically = true

  /**
   * Returns the list of QScripts passed in via -S so that their
   * arguments can be inspected before QScript.script is called.
   * @return Array of QScripts passed in.
   */
  override def getArgumentSources =
    qScriptManager.getValues.asInstanceOf[Array[Class[_]]]

  /**
   * Returns the name of a QScript
   * @return The name of a QScript
   */
  override def getArgumentSourceName(source: Class[_]) =
    qScriptManager.getName(source.asSubclass(classOf[QScript]))

  /**
   * Returns a ScalaCompoundArgumentTypeDescriptor that can parse argument sources into scala collections.
   * @return a ScalaCompoundArgumentTypeDescriptor
   */
  override def getArgumentTypeDescriptors =
    Arrays.asList(new ScalaCompoundArgumentTypeDescriptor)

  /**
   * Loads the QScripts passed in and returns a new QScriptManager than can be used to create them.
   */
  private lazy val qScriptManager = {
    QScriptManager.loadScripts(scripts)
    new QScriptManager
  }
}

/**
 * Entry point of Queue.  Compiles and runs QScripts passed in to the command line.
 */
object QCommandLine {
  /**
   * Main.
   * @param argv Arguments.
   */
  def main(argv: Array[String]) {
    try {
      CommandLineProgram.start(new QCommandLine, argv);
      if (CommandLineProgram.result != 0)
        System.exit(CommandLineProgram.result);
    } catch {
      case e: Exception => CommandLineProgram.exitSystemWithError(e)
    }
  }
}

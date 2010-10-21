package org.broadinstitute.sting.queue

import java.io.File
import java.util.Arrays
import org.broadinstitute.sting.queue.engine.QGraph
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.util._

/**
 * Entry point of Queue.  Compiles and runs QScripts passed in to the command line.
 */
class QCommandLine extends CommandLineProgram with Logging {
  @Input(fullName="script", shortName="S", doc="QScript scala file", required=true)
  @ClassType(classOf[File])
  private var scripts = List.empty[File]

  @Argument(fullName="bsub_all_jobs", shortName="bsub", doc="Use bsub to submit jobs", required=false)
  private var bsubAllJobs = false

  @Argument(fullName="run_scripts", shortName="run", doc="Run QScripts.  Without this flag set only performs a dry run.", required=false)
  private var run = false

  @Argument(fullName="dot_graph", shortName="dot", doc="Outputs the queue graph to a .dot file.  See: http://en.wikipedia.org/wiki/DOT_language", required=false)
  private var dotFile: File = _

  @Argument(fullName="expanded_dot_graph", shortName="expandedDot", doc="Outputs the queue graph of scatter gather to a .dot file.  Otherwise overwrites the dot_graph", required=false)
  private var expandedDotFile: File = _

  @Argument(fullName="start_from_scratch", shortName="startFromScratch", doc="Runs all command line functions even if the outputs were previously output successfully.", required=false)
  private var startFromScratch = false

  @Argument(fullName="for_reals", shortName="forReals", doc="Run QScripts", required=false) @Hidden
  private var runScripts = false

  @Argument(fullName="status",shortName="status",doc="Get status of jobs for the qscript",required=false)
  private var getStatus = false

  @ArgumentCollection
  private val qSettings = new QSettings

  @Argument(fullName="status_email_from", shortName="statusFrom", doc="Email address to send emails from upon completion or on error.", required=false)
  private var statusEmailFrom: String = System.getProperty("user.name") + "@" + SystemUtils.domainName

  @Argument(fullName="status_email_to", shortName="statusTo", doc="Email address to send emails to upon completion or on error.", required=false)
  private var statusEmailTo: List[String] = Nil

  @Argument(fullName="delete_intermediate_outputs", shortName="deleteIntermediates", doc="After a successful run delete the outputs of any Function marked as intermediate.", required=false)
  private var deleteIntermediates = false

  /**
   * Takes the QScripts passed in, runs their script() methods, retrieves their generated
   * functions, and then builds and runs a QGraph based on the dependencies.
   */
  def execute = {

    val qGraph = new QGraph
    qGraph.dryRun = !(run || runScripts)
    qGraph.bsubAllJobs = bsubAllJobs
    qGraph.startFromScratch = startFromScratch
    qGraph.dotFile = dotFile
    qGraph.expandedDotFile = expandedDotFile
    qGraph.qSettings = qSettings
    qGraph.debugMode = debugMode == true
    qGraph.statusEmailFrom = statusEmailFrom
    qGraph.statusEmailTo = statusEmailTo
    qGraph.deleteIntermediates = deleteIntermediates

    val scripts = qScriptManager.createScripts()
    for (script <- scripts) {
      logger.info("Scripting " + qScriptManager.getName(script.getClass.asSubclass(classOf[QScript])))
      loadArgumentsIntoObject(script)
      script.script
      script.functions.foreach(qGraph.add(_))
      logger.info("Added " + script.functions.size + " functions")
    }

    Runtime.getRuntime.addShutdownHook(new Thread {
      /** Cleanup as the JVM shuts down. */
      override def run = {
        qGraph.shutdown()
        ProcessController.shutdown()
        QScriptManager.deleteOutdir()
      }
    })

    if (getStatus)
      qGraph.checkStatus
    else
      qGraph.run

    if (qGraph.hasFailed) {
      logger.info("Done with errors")
      qGraph.logFailed
      1
    } else {
      logger.info("Done")
      0
    }
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

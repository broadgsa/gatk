package org.broadinstitute.sting.queue

import function.QFunction
import java.io.File
import java.util.Arrays
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.engine.{QGraphSettings, QGraph}

/**
 * Entry point of Queue.  Compiles and runs QScripts passed in to the command line.
 */
class QCommandLine extends CommandLineProgram with Logging {
  @Input(fullName="script", shortName="S", doc="QScript scala file", required=true)
  @ClassType(classOf[File])
  private var scripts = List.empty[File]

  @ArgumentCollection
  private val settings = new QGraphSettings

  QFunction.parsingEngine = new ParsingEngine(this)    

  /**
   * Takes the QScripts passed in, runs their script() methods, retrieves their generated
   * functions, and then builds and runs a QGraph based on the dependencies.
   */
  def execute = {

    val qGraph = QCommandLine.qGraph
    qGraph.settings = settings
    qGraph.debugMode = debugMode == true

    val scripts = qScriptManager.createScripts()
    for (script <- scripts) {
      logger.info("Scripting " + qScriptManager.getName(script.getClass.asSubclass(classOf[QScript])))
      loadArgumentsIntoObject(script)
      script.script
      script.functions.foreach(qGraph.add(_))
      logger.info("Added " + script.functions.size + " functions")
    }

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
  private val qGraph = new QGraph


  /**
   * Main.
   * @param argv Arguments.
   */
  def main(argv: Array[String]) {
    Runtime.getRuntime.addShutdownHook(new Thread {
      /** Cleanup as the JVM shuts down. */
      override def run = {
        qGraph.shutdown()
        ProcessController.shutdown()
        QScriptManager.deleteOutdir()
      }
    })

    try {
      CommandLineProgram.start(new QCommandLine, argv);
      if (CommandLineProgram.result != 0)
        System.exit(CommandLineProgram.result);
    } catch {
      case e: Exception => CommandLineProgram.exitSystemWithError(e)
    }
  }
}

package org.broadinstitute.sting.queue

import function.QFunction
import java.io.File
import java.util.Arrays
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.engine.{QGraphSettings, QGraph}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.classloader.PluginManager

/**
 * Entry point of Queue.  Compiles and runs QScripts passed in to the command line.
 */
class QCommandLine extends CommandLineProgram with Logging {
  @Input(fullName="script", shortName="S", doc="QScript scala file", required=true)
  @ClassType(classOf[File])
  private var scripts = List.empty[File]

  @ArgumentCollection
  private val settings = new QGraphSettings

  private val qScriptManager = new QScriptManager
  private val qGraph = new QGraph
  private var qScriptClasses: File = _

  private lazy val pluginManager = {
    qScriptClasses = IOUtils.tempDir("Q-Classes", "", settings.qSettings.tempDirectory)
    qScriptManager.loadScripts(scripts, qScriptClasses)
    new PluginManager[QScript](classOf[QScript], List(qScriptClasses.toURI.toURL))
  }

  QFunction.parsingEngine = new ParsingEngine(this)

  /**
   * Takes the QScripts passed in, runs their script() methods, retrieves their generated
   * functions, and then builds and runs a QGraph based on the dependencies.
   */
  def execute = {
    qGraph.settings = settings
    qGraph.debugMode = debugMode == true

    for (script <- pluginManager.createAllTypes()) {
      logger.info("Scripting " + pluginManager.getName(script.getClass.asSubclass(classOf[QScript])))
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
    pluginManager.getPlugins.toIterable.toArray.asInstanceOf[Array[Class[_]]]

  /**
   * Returns the name of a QScript
   * @return The name of a QScript
   */
  override def getArgumentSourceName(source: Class[_]) =
    pluginManager.getName(source.asSubclass(classOf[QScript]))

  /**
   * Returns a ScalaCompoundArgumentTypeDescriptor that can parse argument sources into scala collections.
   * @return a ScalaCompoundArgumentTypeDescriptor
   */
  override def getArgumentTypeDescriptors =
    Arrays.asList(new ScalaCompoundArgumentTypeDescriptor)

  def shutdown() = {
    qGraph.shutdown()
    if (qScriptClasses != null) IOUtils.tryDelete(qScriptClasses)
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
    val qCommandLine = new QCommandLine

    Runtime.getRuntime.addShutdownHook(new Thread {
      /** Cleanup as the JVM shuts down. */
      override def run = {
        ProcessController.shutdown()
        qCommandLine.shutdown()
      }
    })

    try {
      CommandLineProgram.start(qCommandLine, argv);
      if (CommandLineProgram.result != 0)
        System.exit(CommandLineProgram.result);
    } catch {
      case e: Exception => CommandLineProgram.exitSystemWithError(e)
    } finally {

    }
  }
}

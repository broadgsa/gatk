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

package org.broadinstitute.sting.queue

import function.QFunction
import java.io.File
import java.util.Arrays
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.engine.{QGraphSettings, QGraph}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.classloader.PluginManager
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.utils.io.IOUtils

/**
 * Entry point of Queue.  Compiles and runs QScripts passed in to the command line.
 */
object QCommandLine extends Logging {
  /**
   * Main.
   * @param argv Arguments.
   */
  def main(argv: Array[String]) {
    val qCommandLine = new QCommandLine

    val shutdownHook = new Thread {
      override def run() {
        logger.info("Shutting down jobs. Please wait...")
        qCommandLine.shutdown()
      }
    }

    Runtime.getRuntime.addShutdownHook(shutdownHook)

    try {
      CommandLineProgram.start(qCommandLine, argv)
      try {
        Runtime.getRuntime.removeShutdownHook(shutdownHook)
      } catch {
        case _ => /* ignore, example 'java.lang.IllegalStateException: Shutdown in progress' */
      }
      if (CommandLineProgram.result != 0)
        System.exit(CommandLineProgram.result);
    } catch {
      case e: Exception => CommandLineProgram.exitSystemWithError(e)
    }
  }
}

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
  private var shuttingDown = false

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

    val allQScripts = pluginManager.createAllTypes();
    for (script <- allQScripts) {
      logger.info("Scripting " + pluginManager.getName(script.getClass.asSubclass(classOf[QScript])))
      loadArgumentsIntoObject(script)
      try {
        script.script()
      } catch {
        case e: Exception =>
          throw new UserException.CannotExecuteQScript(script.getClass.getSimpleName + ".script() threw the following exception: " + e, e)
      }
      script.functions.foreach(qGraph.add(_))
      logger.info("Added " + script.functions.size + " functions")
    }

    // Execute the job graph
    qGraph.run()

    // walk over each script, calling onExecutionDone
    for (script <- allQScripts) {
      script.onExecutionDone(qGraph.getFunctionsAndStatus(script.functions), qGraph.success)
      if ( ! settings.disableJobReport ) {
        val jobStringName = (QScriptUtils.?(settings.jobReportFile)).getOrElse(settings.qSettings.jobNamePrefix + ".jobreport.txt")

        if (!shuttingDown) {
          val reportFile = new File(jobStringName)
          logger.info("Writing JobLogging GATKReport to file " + reportFile)
          QJobReport.printReport(qGraph.getFunctionsAndStatus(script.functions), reportFile)

          val pdfFile = new File(jobStringName + ".pdf")
          logger.info("Plotting JobLogging GATKReport to file " + pdfFile)
          QJobReport.plotReport(reportFile, pdfFile)
        }
      }
    }

    if (!qGraph.success) {
      logger.info("Done with errors")
      qGraph.logFailed()
      1
    } else {
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
    shuttingDown = true
    qGraph.shutdown()
    if (qScriptClasses != null) IOUtils.tryDelete(qScriptClasses)
  }
}

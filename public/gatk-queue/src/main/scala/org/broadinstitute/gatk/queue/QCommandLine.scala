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

package org.broadinstitute.gatk.queue

import java.io.File
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.util._
import org.broadinstitute.gatk.queue.engine.{QStatusMessenger, QGraphSettings, QGraph}
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.classloader.PluginManager
import org.broadinstitute.gatk.utils.exceptions.UserException
import org.broadinstitute.gatk.utils.io.IOUtils
import org.broadinstitute.gatk.utils.help.ApplicationDetails
import java.util.Arrays
import org.apache.commons.io.FilenameUtils

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
        qCommandLine.shutdown()
      } catch {
        case e: Exception => /* ignore, example 'java.lang.IllegalStateException: Shutdown in progress' */
      }
      if (CommandLineProgram.result != 0)
        System.exit(CommandLineProgram.result)
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
  var scripts = Seq.empty[File]

  @ArgumentCollection
  val settings = new QGraphSettings

  private val qScriptManager = new QScriptManager
  private val qGraph = new QGraph
  private var qScriptClasses: File = _
  private var shuttingDown = false

  private lazy val qScriptPluginManager = {
    qScriptClasses = IOUtils.tempDir("Q-Classes-", "", settings.qSettings.tempDirectory)
    qScriptManager.loadScripts(scripts, qScriptClasses)
    new PluginManager[QScript](qPluginType, Seq(qScriptClasses.toURI.toURL))
  }

  private lazy val qCommandPlugin = {
    new PluginManager[QCommandPlugin](classOf[QCommandPlugin])
  }

  private lazy val allCommandPlugins = qCommandPlugin.createAllTypes()

  private lazy val qPluginType: Class[_ <: QScript] = {
    allCommandPlugins.map(_.qScriptClass).headOption.getOrElse(classOf[QScript])
  }

  /**
   * Takes the QScripts passed in, runs their script() methods, retrieves their generated
   * functions, and then builds and runs a QGraph based on the dependencies.
   */
  def execute = {
    var success = false
    var result = 1
    var functionsAndStatusSize = 0
    try {
      ClassFieldCache.parsingEngine = this.parser

      if (settings.qSettings.runName == null)
        settings.qSettings.runName = FilenameUtils.removeExtension(scripts.head.getName)
      if (IOUtils.isDefaultTempDir(settings.qSettings.tempDirectory))
        settings.qSettings.tempDirectory = IOUtils.absolute(settings.qSettings.runDirectory, ".queue/tmp")
      qGraph.initializeWithSettings(settings)

      for (commandPlugin <- allCommandPlugins) {
        loadArgumentsIntoObject(commandPlugin)
      }

      for (commandPlugin <- allCommandPlugins) {
        if (commandPlugin.statusMessenger != null)
          commandPlugin.statusMessenger.started()
      }

      qGraph.messengers = allCommandPlugins.filter(_.statusMessenger != null).map(_.statusMessenger).toSeq

      // TODO: Default command plugin argument?
      val remoteFileConverter = (
        for (commandPlugin <- allCommandPlugins if (commandPlugin.remoteFileConverter != null))
        yield commandPlugin.remoteFileConverter
        ).headOption.getOrElse(null)

      if (remoteFileConverter != null)
        loadArgumentsIntoObject(remoteFileConverter)

      val allQScripts = qScriptPluginManager.createAllTypes()
      for (script <- allQScripts) {
        logger.info("Scripting " + qScriptPluginManager.getName(script.getClass.asSubclass(classOf[QScript])))
        loadArgumentsIntoObject(script)
        allCommandPlugins.foreach(_.initScript(script))
        // TODO: Pulling inputs can be time/io expensive! Some scripts are using the files to generate functions-- even for dry runs-- so pull it all down for now.
        //if (settings.run)
        script.pullInputs()
        script.qSettings = settings.qSettings
        try {
          script.script()
        } catch {
          case e: Exception =>
            throw new UserException.CannotExecuteQScript(script.getClass.getSimpleName + ".script() threw the following exception: " + e, e)
        }

        if (remoteFileConverter != null) {
          if (remoteFileConverter.convertToRemoteEnabled)
            script.mkRemoteOutputs(remoteFileConverter)
        }

        script.functions.foreach(qGraph.add(_))
        logger.info("Added " + script.functions.size + " functions")
      }
      // Execute the job graph
      qGraph.run()

      val functionsAndStatus = qGraph.getFunctionsAndStatus

      // walk over each script, calling onExecutionDone
      for (script <- allQScripts) {
        val scriptFunctions = functionsAndStatus.filterKeys(f => script.functions.contains(f))
        script.onExecutionDone(scriptFunctions, success)
      }
      functionsAndStatusSize = functionsAndStatus.size

      // write the final complete job report
      logger.info("Writing final jobs report...")
      qGraph.writeJobsReport()

      if (qGraph.success) {
        if (settings.run) {
          allQScripts.foreach(_.pushOutputs())
          for (commandPlugin <- allCommandPlugins)
            if (commandPlugin.statusMessenger != null) {
              val allInputs = allQScripts.map(_.remoteInputs)
              val allOutputs = allQScripts.map(_.remoteOutputs)
              commandPlugin.statusMessenger.done(allInputs, allOutputs)
            }
        }
        success = true
        result = 0
      }
    } finally {
      if (!success) {
        logger.info("Done with errors")
        qGraph.logFailed()
        if (settings.run) {
          for (commandPlugin <- allCommandPlugins)
            if (commandPlugin.statusMessenger != null)
              commandPlugin.statusMessenger.exit("Done with errors: %s".format(qGraph.formattedStatusCounts))
        }
      }
    }
    if (success)
      logger.info("Script completed successfully with %d total jobs".format(functionsAndStatusSize))
    else
      logger.info("Script failed: %s".format(qGraph.formattedStatusCounts))
    result
  }

  /**
   * Returns true as QScripts are located and compiled.
   * @return true
   */
  override def canAddArgumentsDynamically = true

  /**
   * Returns the list of QScripts passed in via -S and other plugins
   * so that their arguments can be inspected before QScript.script is called.
   * @return Array of dynamic sources
   */
  override def getArgumentSources = {
    var plugins = Seq.empty[Class[_]]
    plugins ++= qScriptPluginManager.getPlugins
    plugins ++= qCommandPlugin.getPlugins
    plugins.toArray
  }

  /**
   * Returns the name of a script/plugin
   * @return The name of a script/plugin
   */
  override def getArgumentSourceName(source: Class[_]) = {
    if (classOf[QScript].isAssignableFrom(source))
      qScriptPluginManager.getName(source.asSubclass(classOf[QScript]))
    else if (classOf[QCommandPlugin].isAssignableFrom(source))
      qCommandPlugin.getName(source.asSubclass(classOf[QCommandPlugin]))
    else
      null
  }

  /**
   * Returns a ScalaCompoundArgumentTypeDescriptor that can parse argument sources into scala collections.
   * @return a ScalaCompoundArgumentTypeDescriptor
   */
  override def getArgumentTypeDescriptors =
    Arrays.asList(new ScalaCompoundArgumentTypeDescriptor)

  override def getApplicationDetails : ApplicationDetails = {
    new ApplicationDetails(createQueueHeader(),
                           Seq.empty[String],
                           ApplicationDetails.createDefaultRunningInstructions(getClass.asInstanceOf[Class[CommandLineProgram]]),
                           "")
  }

  private def createQueueHeader() : Seq[String] = {
    Seq(String.format("Queue v%s, Compiled %s", CommandLineProgram.getVersionNumber, CommandLineProgram.getBuildTime),
         "Copyright (c) 2012 The Broad Institute",
         "For support and documentation go to http://www.broadinstitute.org/gatk")
  }

  def shutdown() = {
    shuttingDown = true
    qGraph.shutdown()
    if (qScriptClasses != null) IOUtils.tryDelete(qScriptClasses)
  }
}

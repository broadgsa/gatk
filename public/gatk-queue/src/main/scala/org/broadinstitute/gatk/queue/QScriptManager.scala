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

import scala.tools.nsc.{Global, Settings}
import scala.tools.nsc.io.PlainFile
import org.broadinstitute.gatk.queue.util.Logging
import collection.JavaConversions._
import java.io.File
import scala.tools.nsc.reporters.AbstractReporter
import java.lang.String
import org.apache.log4j.Level
import org.broadinstitute.gatk.queue.util.TextFormatUtils._
import org.broadinstitute.gatk.utils.classloader.JVMUtils
import scala.reflect.internal.util.{FakePos, NoPosition, Position, StringOps}
import org.broadinstitute.gatk.utils.exceptions.UserException

/**
 * Plugin manager for QScripts which loads QScripts into the current class loader.
 */
class QScriptManager() extends Logging {
  /**
   * Compiles and loads the scripts in the files into the current classloader.
   * Heavily based on scala/src/compiler/scala/tools/ant/Scalac.scala
   */
  def loadScripts(scripts: Seq[File], tempDir: File) {
    // Make sure the scripts actually exist.
    scripts.foreach{
        file => if( !file.exists()) throw new UserException.CouldNotReadInputFile(file, "it does not exist.")
    }

    if (scripts.size > 0) {
      val settings = new Settings((error: String) => logger.error(error))
      settings.deprecation.value = true
      settings.outdir.value = tempDir.getPath

      // Set the classpath to the current class path.
      JVMUtils.getClasspathURLs.foreach(url => {
          settings.bootclasspath.append(url.getPath)
          settings.classpath.append(url.getPath)
      })

      val reporter = new QScriptManager.Log4JReporter(settings)

      val compiler = new Global(settings, reporter)
      val run = new compiler.Run

      logger.info("Compiling %s QScript%s".format(scripts.size, plural(scripts.size)))
      logger.debug("Compilation directory: " + settings.outdir.value)
      run.compileFiles(scripts.toList.map(new PlainFile(_)))

      reporter.printSummary()
      if (reporter.hasErrors) {
        val msg = "Compile of %s failed with %d error%s".format(
          scripts.mkString(", "), reporter.ERROR.count, plural(reporter.ERROR.count))
        throw new QException(msg)
      }
      else if (reporter.WARNING.count > 0)
        logger.warn("Compile succeeded with %d warning%s".format(
          reporter.WARNING.count, plural(reporter.WARNING.count)))
      else
        logger.info("Compilation complete")
    }
  }
}

/**
 * Plugin manager for QScripts which loads QScripts into the current classloader.
 */
object QScriptManager extends Logging {

  /**
   * NSC (New Scala Compiler) reporter which logs to Log4J.
   * Heavily based on scala/src/compiler/scala/tools/nsc/reporters/ConsoleReporter.scala
   */
  private class Log4JReporter(val settings: Settings) extends AbstractReporter {
    def displayPrompt() { throw new UnsupportedOperationException("Unable to prompt the user.  Prompting should be off.") }

    /**
     * Displays the message at position with severity.
     * @param posIn Position of the event in the file that generated the message.
     * @param msg Message to display.
     * @param severity Severity of the event.
     */
    def display(posIn: Position, msg: String, severity: Severity) {
      severity.count += 1
      val level = severity match {
        case INFO => Level.INFO
        case WARNING => Level.WARN
        case ERROR => Level.ERROR
      }
      val pos = if (posIn eq null) NoPosition
                else if (posIn.isDefined) posIn.inUltimateSource(posIn.source)
                else posIn
      pos match {
        case FakePos(fmsg) =>
          printMessage(level, fmsg+" "+msg)
        case NoPosition =>
          printMessage(level, msg)
        case _ =>
          val file = pos.source.file
          printMessage(level, file.name+":"+pos.line+": "+msg)
          printSourceLine(level, pos)
      }
    }

    /**
     * Prints a summary count of warnings and errors.
     */
    def printSummary() {
      if (WARNING.count > 0)
        printMessage(Level.WARN, StringOps.countElementsAsString(WARNING.count, "warning") + " found")
      if (ERROR.count > 0)
        printMessage(Level.ERROR, StringOps.countElementsAsString(ERROR.count, "error") + " found")
    }

    /**
     * Prints the source code line of an event followed by a pointer within the line to the error.
     * @param level Severity level.
     * @param pos Position in the file of the event.
     */
    private def printSourceLine(level: Level, pos: Position) {
      printMessage(level, pos.lineContent.stripLineEnd)
      printColumnMarker(level, pos)
    }

    /**
     * Prints the column marker of the given position.
     * @param level Severity level.
     * @param pos Position in the file of the event.
     */
    private def printColumnMarker(level: Level, pos: Position) {
      if (pos.isDefined) { printMessage(level, " " * (pos.column - 1) + "^") }
    }

    /**
     * Prints the message at the severity level.
     * @param level Severity level.
     * @param message Message content.
     */
    private def printMessage(level: Level, message: String) {
      logger.log(level, message)
    }
  }
}

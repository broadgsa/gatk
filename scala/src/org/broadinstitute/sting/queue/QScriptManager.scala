package org.broadinstitute.sting.queue

import scala.tools.nsc.{Global, Settings}
import scala.tools.nsc.io.PlainFile
import org.broadinstitute.sting.queue.util.{Logging, IOUtils}
import collection.JavaConversions
import java.io.File
import scala.tools.nsc.reporters.AbstractReporter
import java.lang.String
import org.apache.log4j.Level
import scala.tools.nsc.util.{FakePos, NoPosition, Position}
import org.broadinstitute.sting.utils.classloader.{PackageUtils, PluginManager}

/**
 * Plugin manager for QScripts which loads QScripts into the current class loader.
 */
class QScriptManager extends PluginManager[QScript](classOf[QScript], "QScript", "Script") with Logging {

  /**
   * Returns the list of QScripts classes found in the classpath.
   * @return QScripts classes found in the classpath.
   */
  def getValues = {
    if (logger.isDebugEnabled) {
      JavaConversions.asMap(this.pluginsByName)
              .foreach{case (name, clazz) => logger.debug("Found QScript %s: %s".format(name, clazz))}
    }
    JavaConversions.asIterable(this.pluginsByName.values).toArray
  }

  /**
   * Creates the QScripts for all values found in the classpath.
   * @return QScripts found in the classpath.
   */
  def createScripts() = getValues.map(_.newInstance.asInstanceOf[QScript])
}

/**
 * Plugin manager for QScripts which loads QScripts into the current classloader.
 */
object QScriptManager extends Logging {
  /**
   * Compiles and loads the scripts in the files into the current classloader.
   * Heavily based on scala/src/compiler/scala/tools/ant/Scalac.scala
   * @param scripts Scala classes to compile.
   */
  def loadScripts(scripts: List[File]) {
    if (scripts.size > 0) {

      val settings = new Settings((error: String) => logger.error(error))
      val outdir = IOUtils.tempDir("Q-classes")
      settings.outdir.value = outdir.getPath

      // Set the classpath to the current class path.
      JavaConversions.asSet(PackageUtils.getClassPathURLs).foreach(url => settings.classpath.append(url.getPath))

      val reporter = new Log4JReporter(settings)

      val compiler = new Global(settings, reporter)
      val run = new compiler.Run

      logger.info("Compiling %s QScript%s".format(scripts.size, plural(scripts.size)))
      logger.debug("Compilation directory: " + settings.outdir.value)
      run.compileFiles(scripts.map(new PlainFile(_)))

      reporter.printSummary()
      if (reporter.hasErrors) {
        val msg = "Compile failed with %d error%s".format(
          reporter.ERROR.count, plural(reporter.ERROR.count))
        throw new QException(msg)
      }
      else if (reporter.WARNING.count > 0)
        logger.warn("Compile succeeded with %d warning%s".format(
          reporter.WARNING.count, plural(reporter.WARNING.count)))
      else
        logger.info("Compilation complete")

      // Add the new compilation output directory to the classpath.
      PackageUtils.addClasspath(outdir.toURI.toURL)
    }
  }

  /**
   * Returns the string "s" if x is greater than 1.
   * @param x Value to test.
   * @return "s" if x is greater than one else "".
   */
  private def plural(x: Int) = if (x > 1) "s" else ""

  /**
   * NSC (New Scala Compiler) reporter which logs to Log4J.
   * Heavily based on scala/src/compiler/scala/tools/nsc/reporters/ConsoleReporter.scala
   */
  private class Log4JReporter(val settings: Settings) extends AbstractReporter {
    def displayPrompt = throw new UnsupportedOperationException("Unable to prompt the user.  Prompting should be off.")

    /**
     * Displays the message at position with severity.
     * @param posIn Position of the event in the file that generated the message.
     * @param msg Message to display.
     * @param severity Severity of the event.
     */
    def display(posIn: Position, msg: String, severity: Severity) = {
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
          val buf = new StringBuilder(msg)
          val file = pos.source.file
          printMessage(level, file.name+":"+pos.line+": "+msg)
          printSourceLine(level, pos)
      }
    }

    /**
     * Prints a summary count of warnings and errors.
     */
    def printSummary() = {
      if (WARNING.count > 0)
        printMessage(Level.WARN, countElementsAsString(WARNING.count, "warning") + " found")
      if (ERROR.count > 0)
        printMessage(Level.ERROR, countElementsAsString(ERROR.count, "error") + " found")
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
    private def printColumnMarker(level: Level, pos: Position) =
      if (pos.isDefined) { printMessage(level, " " * (pos.column - 1) + "^") }

    /**
     * Prints the message at the severity level.
     * @param level Severity level.
     * @param message Message content.
     */
    private def printMessage(level: Level, message: String) = {
      logger.log(level, message)
    }
  }
}

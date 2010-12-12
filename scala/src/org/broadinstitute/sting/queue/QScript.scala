package org.broadinstitute.sting.queue

import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.QFunction

/**
 * Defines a Queue pipeline as a collection of CommandLineFunctions.
 */
trait QScript extends Logging {
  // Type aliases so users don't have to import
  type File = java.io.File
  type Input = org.broadinstitute.sting.commandline.Input
  type Output = org.broadinstitute.sting.commandline.Output
  type Argument = org.broadinstitute.sting.commandline.Argument
  type ArgumentCollection = org.broadinstitute.sting.commandline.ArgumentCollection
  type CommandLineFunction = org.broadinstitute.sting.queue.function.CommandLineFunction
  type InProcessFunction = org.broadinstitute.sting.queue.function.InProcessFunction
  type ScatterGatherableFunction = org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
  type Gather = org.broadinstitute.sting.queue.function.scattergather.Gather
  type SimpleTextGatherFunction = org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction

  /**
   * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
   */
  def script: Unit

  /**
   * The command line functions that will be executed for this QScript.
   */
  var functions = List.empty[QFunction]

  /**
   * Exchanges the extension on a file.
   * @param file File to look for the extension.
   * @param oldExtension Old extension to strip off, if present.
   * @param newExtension New extension to append.
   * @return new File with the new extension in the current directory.
   */
  protected def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getName.stripSuffix(oldExtension) + newExtension)

  /**
   * Exchanges the extension on a file.
   * @param dir New directory for the file.
   * @param file File to look for the extension.
   * @param oldExtension Old extension to strip off, if present.
   * @param newExtension New extension to append.
   * @return new File with the new extension in dir.
   */
  protected def swapExt(dir: String, file: File, oldExtension: String, newExtension: String) =
    new File(dir, file.getName.stripSuffix(oldExtension) + newExtension)

  /**
   * Exchanges the extension on a file.
   * @param dir New directory for the file.
   * @param file File to look for the extension.
   * @param oldExtension Old extension to strip off, if present.
   * @param newExtension New extension to append.
   * @return new File with the new extension in dir.
   */
  protected def swapExt(dir: File, file: File, oldExtension: String, newExtension: String) =
    new File(dir, file.getName.stripSuffix(oldExtension) + newExtension)

  /**
   * Adds one or more command line functions to be run.
   * @param functions Functions to add.
   */
  def add(functions: QFunction*) = {
    functions.foreach(function => function.addOrder = QScript.nextAddOrder)
    this.functions ++= functions
  }

  def addAll(functions: List[QFunction] ) = {
    functions.foreach( f => add(f) )
  }
}

object QScript {
  private var addOrder = 0
  private def nextAddOrder = {
    addOrder += 1
    List(addOrder)
  }
}

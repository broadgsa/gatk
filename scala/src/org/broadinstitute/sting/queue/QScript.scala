package org.broadinstitute.sting.queue

import org.broadinstitute.sting.queue.util.Logging

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
  type ScatterGatherableFunction = org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
  type Scatter = org.broadinstitute.sting.queue.function.scattergather.Scatter
  type Gather = org.broadinstitute.sting.queue.function.scattergather.Gather
  type SimpleTextGatherFunction = org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction

  /**
   * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
   */
  def script: Unit

  /**
   * The command line functions that will be executed for this QScript.
   */
  var functions = List.empty[CommandLineFunction]

  /**
   * Exchanges the extension on a file.
   */
  protected def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getName.stripSuffix(oldExtension) + newExtension)

  /**
   * Adds one or more command line functions to be run.
   */
  def add(functions: CommandLineFunction*) = this.functions ++= List(functions:_*)

}

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

import engine.JobRunInfo
import org.broadinstitute.sting.queue.function.QFunction
import annotation.target.field
import io.Source
import util.{StringFileConversions, PrimitiveOptionConversions, Logging}

/**
 * Defines a Queue pipeline as a collection of CommandLineFunctions.
 */
trait QScript extends Logging with PrimitiveOptionConversions with StringFileConversions  {

  // Type aliases so users don't have to import
  type File = java.io.File
  type CommandLineFunction = org.broadinstitute.sting.queue.function.CommandLineFunction
  type InProcessFunction = org.broadinstitute.sting.queue.function.InProcessFunction
  type ScatterGatherableFunction = org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
  type SimpleTextGatherFunction = org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction

  // Make sure annotations can be used in class constructors but target the fields
  // ex: class MyClass(@Input var myVar: File) {}
  // This was implicitly enabled in 2.8.0-RC2 and then updated to this new syntax:
  // http://lampsvn.epfl.ch/trac/scala/ticket/3596
  // http://lampsvn.epfl.ch/trac/scala/ticket/3421
  type Input = org.broadinstitute.sting.commandline.Input @field
  type Output = org.broadinstitute.sting.commandline.Output @field
  type Argument = org.broadinstitute.sting.commandline.Argument @field
  type ArgumentCollection = org.broadinstitute.sting.commandline.ArgumentCollection @field
  type Gather = org.broadinstitute.sting.commandline.Gather @field

  /**
   * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
   */
  def script()

  /**
   * A default handler for the onExecutionDone() function.  By default this doesn't do anything
   * except print out a fine status message.
   */
  def onExecutionDone(jobs: Map[QFunction, JobRunInfo], success: Boolean) {
    logger.info("Script %s with %d total jobs".format(if (success) "completed successfully" else "failed", jobs.size))
    // this is too much output
    // for ( (f, info) <- jobs ) logger.info("  %s %s".format(f.jobName, info))
  }

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

  def addAll(functions: List[QFunction]) {
    functions.foreach( f => add(f) )
  }

  def extractFileEntries(in: File): List[File] = Source.fromFile(in).getLines().toList
}

object QScript {
  private var addOrder = 0
  private def nextAddOrder = {
    addOrder += 1
    List(addOrder)
  }
}

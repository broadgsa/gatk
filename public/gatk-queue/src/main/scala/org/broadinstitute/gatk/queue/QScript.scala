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

import engine.JobRunInfo
import org.broadinstitute.gatk.queue.function.QFunction
import annotation.target.field
import util._
import org.broadinstitute.gatk.utils.commandline.ArgumentSource

/**
 * Defines a Queue pipeline as a collection of CommandLineFunctions.
 */
trait QScript extends Logging with PrimitiveOptionConversions with StringFileConversions  {

  // Type aliases so users don't have to import
  type File = java.io.File
  type CommandLineFunction = org.broadinstitute.gatk.queue.function.CommandLineFunction
  type InProcessFunction = org.broadinstitute.gatk.queue.function.InProcessFunction
  type ScatterGatherableFunction = org.broadinstitute.gatk.queue.function.scattergather.ScatterGatherableFunction
  type SimpleTextGatherFunction = org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction

  // Make sure annotations can be used in class constructors but target the fields
  // ex: class MyClass(@Input var myVar: File) {}
  // This was implicitly enabled in 2.8.0-RC2 and then updated to this new syntax:
  // http://lampsvn.epfl.ch/trac/scala/ticket/3596
  // http://lampsvn.epfl.ch/trac/scala/ticket/3421
  type Input = org.broadinstitute.gatk.utils.commandline.Input @field
  type Output = org.broadinstitute.gatk.utils.commandline.Output @field
  type Argument = org.broadinstitute.gatk.utils.commandline.Argument @field
  type ArgumentCollection = org.broadinstitute.gatk.utils.commandline.ArgumentCollection @field
  type Gather = org.broadinstitute.gatk.utils.commandline.Gather @field

  /**
   * Default settings for QFunctions
   */
  var qSettings: QSettings = _

  /**
   * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
   */
  def script()

  /**
   * A default handler for the onExecutionDone() function.  By default this doesn't do anything
   */
  def onExecutionDone(jobs: Map[QFunction, JobRunInfo], success: Boolean) {
  }

  /**
   * The command line functions that will be executed for this QScript.
   */
  var functions = Seq.empty[QFunction]

  /**
   * Exchanges the extension on a file.
   * @param file File to look for the extension.
   * @param oldExtension Old extension to strip off, if present.
   * @param newExtension New extension to append.
   * @return new File with the new extension in the current directory.
   */
  protected def swapExt(file: File, oldExtension: String, newExtension: String) = QScriptUtils.swapExt(file, oldExtension, newExtension)

  /**
   * Exchanges the extension on a file.
   * @param dir New directory for the file.
   * @param file File to look for the extension.
   * @param oldExtension Old extension to strip off, if present.
   * @param newExtension New extension to append.
   * @return new File with the new extension in dir.
   */
  protected def swapExt(dir: File, file: File, oldExtension: String, newExtension: String) = QScriptUtils.swapExt(dir, file, oldExtension, newExtension)

  /**
   * Adds one or more command line functions to be run.
   * @param functions Functions to add.
   */
  def add(functions: QFunction*) {
    functions.foreach(function => function.addOrder = QScript.nextAddOrder)
    this.functions ++= functions
  }

  def addAll(functions: Traversable[QFunction]) {
    functions.foreach( f => add(f) )
  }

  /**
   * Convert all @Output files to remote output files.
   * @param remoteFileConverter Converter for files to remote files.
   */
  def mkRemoteOutputs(remoteFileConverter: RemoteFileConverter) {
    for (field <- outputFields) {
      val fieldFile = ClassFieldCache.getFieldFile(this, field)
      if (fieldFile != null && !fieldFile.isInstanceOf[RemoteFile]) {
        val fieldName = ClassFieldCache.fullName(field)
        val remoteFile = remoteFileConverter.convertToRemote(fieldFile, fieldName)
        ClassFieldCache.setFieldValue(this, field, remoteFile)
      }
    }
  }

  /**
   * Pull all remote files to the local disk
   */
  def pullInputs() {
    val inputs = ClassFieldCache.getFieldFiles(this, inputFields)
    for (remoteFile <- filterRemoteFiles(inputs)) {
      logger.info("Pulling %s from %s".format(remoteFile.getAbsolutePath, remoteFile.remoteDescription))
      remoteFile.pullToLocal()
    }
  }

  /**
   * Push all remote files from the local disk
   */
  def pushOutputs() {
    val outputs = ClassFieldCache.getFieldFiles(this, outputFields)
    for (remoteFile <- filterRemoteFiles(outputs)) {
      logger.info("Pushing %s to %s".format(remoteFile.getAbsolutePath, remoteFile.remoteDescription))
      remoteFile.pushToRemote()
    }
  }

  private def filterRemoteFiles(fields: Seq[File]): Seq[RemoteFile] =
    fields.filter(field => field != null && field.isInstanceOf[RemoteFile]).map(_.asInstanceOf[RemoteFile])
  /**
   * @return the inputs or null if there are no inputs
   */
  def remoteInputs: AnyRef = null

  /**
   * @return the outputs or null if there are no outputs
   */
  def remoteOutputs: AnyRef = null

  /** The complete list of fields. */
  def functionFields: Seq[ArgumentSource] = ClassFieldCache.classFunctionFields(this.getClass)
  /** The @Input fields. */
  def inputFields: Seq[ArgumentSource] = ClassFieldCache.classInputFields(this.getClass)
  /** The @Output fields. */
  def outputFields: Seq[ArgumentSource] = ClassFieldCache.classOutputFields(this.getClass)
  /** The @Argument fields. */
  def argumentFields: Seq[ArgumentSource] = ClassFieldCache.classArgumentFields(this.getClass)
}

object QScript {
  private var addOrder = 0
  private def nextAddOrder = {
    addOrder += 1
    Seq(addOrder)
  }

  /**
   * Resets the add order back to zero. Useful for testing purposes.
   */
  def resetAddOrder() {
    addOrder = 0
  }

}

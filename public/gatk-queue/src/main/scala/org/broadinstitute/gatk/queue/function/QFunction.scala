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

package org.broadinstitute.gatk.queue.function

import java.io.File
import java.lang.annotation.Annotation
import org.broadinstitute.gatk.utils.commandline._
import org.broadinstitute.gatk.queue.{QException, QSettings}
import java.lang.IllegalStateException
import org.broadinstitute.gatk.queue.util._
import org.broadinstitute.gatk.utils.io.IOUtils
import scala.language.reflectiveCalls

/**
 * The base interface for all functions in Queue.
 * Inputs and outputs are specified as Sets of values.
 * Inputs are matched to other outputs by using .equals()
 */
trait QFunction extends Logging with QJobReport {
  /**
   * A short description of what this class of function does.
   * By default does not include the output specific to this function.
   * See shortDescription for a description of what this instance of the function outputs.
   */
  var analysisName: String = "<function>"

  /**
   * The name name of the job, must be file system safe and unique to the graph.
   * Defaults to "runName-<order_function_added>".
   * Use shortDescription for an alternative that is display friendly.
   */
  var jobName: String = _

  /** Default settings */
  var qSettings: QSettings = _

  /** Directory to run the command in. */
  var commandDirectory: File = new File(".")

  /** Temporary directory to write any files. Must be network accessible. */
  var jobTempDir: File = null

  /**
   * Local path available on all machines to store LOCAL temporary files. Not an @Input,
   * nor an @Output. Currently only used for local intermediate files for composite jobs.
   * Needs to be an annotated field so that it's mutated during cloning.
   */
  @Argument(doc="Local path available on all machines to store LOCAL temporary files.")
  var jobLocalDir: File = _

  /** Order the function was added to the graph. */
  var addOrder: Seq[Int] = Nil

  /** Job priority */
  var jobPriority: Option[Int] = None

  /** Whether a job is restartable */
  var jobRestartable = true

  /**
   * A callback for modifying the run.
   * NOTE: This function is for ADVANCED use only and is unsupported.
   */
  var updateJobRun: PartialFunction[Any,Unit] = null

  /**
   * If true, unless another unfinished function is dependent on this function,
   * this function will NOT be run even if the outputs have not been created.
   */
  var isIntermediate = false

  // -------------------------------------------------------
  //
  // job run information
  //
  // -------------------------------------------------------

  /**
   * Copies settings from this function to another function.
   * @param function QFunction to copy values to.
   */
  override def copySettingsTo(function: QFunction) {
    function.qSettings = this.qSettings
    function.commandDirectory = this.commandDirectory
    function.jobTempDir = this.jobTempDir
    function.jobLocalDir = this.jobLocalDir
    function.addOrder = this.addOrder
    function.jobPriority = this.jobPriority
    function.jobRestartable = this.jobRestartable
    function.updateJobRun = this.updateJobRun
    function.isIntermediate = this.isIntermediate
    function.reportGroup = this.reportGroup
    function.reportFeatures = this.reportFeatures
  }

  /** File to redirect any output.  Defaults to <jobName>.out */
  var jobOutputFile: File = _

  /** File to redirect any errors.  Defaults to <jobName>.out */
  var jobErrorFile: File = _

  /** Errors (if any) from the last failed run of jobErrorFiles. */
  @Argument(doc="Job error lines", required=false)
  var jobErrorLines: Seq[String] = Nil

  /**
   * The number of times this function has previously been run.
   */
  @Argument(doc="Job retries", required=false)
  var retries = 0

  /** Change settings for the next run. Retries will be set to the number of times the function was run and jobErrorLines may contain the error text. */
  def setupRetry() {
  }

  /**
   * Description of this command line function.
   */
  def description: String = "%s: %s > %s".format(analysisName, inputs, outputs)

  /**
   * A short description of the function.
   */
  def shortDescription = {
    firstOutput match {
      case file: File => analysisName + ": " + file.getName
      case _ => analysisName
    }
  }
  
  /**
   * The name of the job as submitted to the job runner
   */
  def jobRunnerJobName = shortDescription

  /**
   * Returns true if the function is done.
   */
  def isDone: Boolean = {
    val files = doneOutputs
    if (files.size == 0)
      throw new IllegalStateException("Function should have at least one output: " + analysisName)
    files.forall(_.exists)
  }

  /**
   * Returns true if the function has failed.
   */
  def isFail: Boolean = {
    val files = failOutputs
    if (files.size == 0)
      throw new IllegalStateException("Function should have at least one output: " + analysisName)
    files.exists(_.exists)
  }

  /**
   * Returns files to track for hidden done/fail files.
   * @return Seq[String] files.
   */
  protected def statusPaths = {
    var paths = outputs
    paths :+= jobOutputFile
    if (jobErrorFile != null)
      paths :+= jobErrorFile
    paths
  }

  /**
   * Returns prefixes for hidden done/fail files.
   * @return prefixes.
   */
  private def statusPrefixes = statusPaths.
    filter(file => !IOUtils.isSpecialFile(file)).
    map(file => file.getParentFile + "/." + file.getName)

  /**
   * Returns the output files for this function.
   * @return outputs for this function.
   */
  def doneOutputs: Seq[File] = statusPrefixes.map(path => new File(path + ".done"))

  /**
   * Returns the output files for this function.
   * @return outputs for this function.
   */
  def failOutputs: Seq[File] = statusPrefixes.map(path => new File(path + ".fail"))

  /** The complete list of fields on this CommandLineFunction. */
  def functionFields: Seq[ArgumentSource] = ClassFieldCache.classFunctionFields(this.functionFieldClass)
  /** The @Input fields on this CommandLineFunction. */
  def inputFields: Seq[ArgumentSource] = ClassFieldCache.classInputFields(this.functionFieldClass)
  /** The @Output fields on this CommandLineFunction. */
  def outputFields: Seq[ArgumentSource] = ClassFieldCache.classOutputFields(this.functionFieldClass)
  /** The @Argument fields on this CommandLineFunction. */
  def argumentFields: Seq[ArgumentSource] = ClassFieldCache.classArgumentFields(this.functionFieldClass)

  /**
   * Returns the class that should be used for looking up fields.
   */
  protected def functionFieldClass = this.getClass

  /**
   * Returns the input files for this function.
   * @return inputs for this function.
   */
  def inputs: Seq[File] = getFieldFiles(inputFields)

  /**
   * Returns the output files for this function.
   * @return outputs for this function.
   */
  def outputs: Seq[File] = getFieldFiles(outputFields)

  /**
   * Returns the first output file.
   * @return first output for this function.
   */
  def firstOutput: File = outputs.headOption.getOrElse(null)

  /**
   * Returns the set of directories where files may be written.
   */
  def outputDirectories = {
    var dirs = Set.empty[File]
    dirs += commandDirectory
    dirs += jobTempDir
    dirs += jobLocalDir
    dirs += jobOutputFile.getParentFile
    if (jobErrorFile != null)
      dirs += jobErrorFile.getParentFile
    dirs ++= outputs.map(_.getParentFile)
    dirs
  }

  /**
   * Deletes the log files for this function.
   */
  def deleteLogs() = {
    IOUtils.tryDelete(jobOutputFile)
    if (jobErrorFile != null)
      IOUtils.tryDelete(jobErrorFile)
  }

  /**
   * Deletes the output files and all the status files for this function.
   */
  def deleteOutputs() {
    outputs.filter(file => !IOUtils.isSpecialFile(file)).foreach(file => IOUtils.tryDelete(file))
    doneOutputs.foreach(file => IOUtils.tryDelete(file))
    failOutputs.foreach(file => IOUtils.tryDelete(file))
  }

  /**
   * Creates the output directories for this function if it doesn't exist.
   */
  def mkOutputDirectories() {
    outputDirectories.foreach(dir => {
      if (!dir.exists && !dir.mkdirs)
        throw new QException("Unable to create directory: " + dir)
    })
  }

  /**
   * Returns fields that do not have values which are required.
   * @return Seq[String] names of fields missing values.
   */
  def missingFields: Seq[String] = {
    val missingInputs = missingFields(inputFields, classOf[Input])
    val missingOutputs = missingFields(outputFields, classOf[Output])
    val missingArguments = missingFields(argumentFields, classOf[Argument])
    (missingInputs ++ missingOutputs ++ missingArguments).distinct.sorted
  }

  /**
   * Returns fields that do not have values which are required.
   * @param sources Fields to check.
   * @param annotation Annotation.
   * @return names of fields missing values.
   */
  private def missingFields(sources: Seq[ArgumentSource], annotation: Class[_ <: Annotation]): Seq[String] = {
    var missing: Seq[String] = Nil
    for (source <- sources) {
      if (isRequired(source, annotation))
        if (!hasFieldValue(source))
          if (!exclusiveOf(source, annotation).exists(otherSource => hasFieldValue(otherSource)))
            missing :+= "@%s: %s - %s".format(annotation.getSimpleName, source.field.getName, doc(source, annotation))
    }
    missing
  }

  /**
   *  Gets the files from the fields.  The fields must be a File, a FileExtension, or a Seq or Set of either.
   * @param fields Fields to get files.
   * @return for the fields.
   */
  private def getFieldFiles(fields: Seq[ArgumentSource]): Seq[File] = {
    var files: Seq[File] = Nil
    for (field <- fields)
      files ++= getFieldFiles(field)
    files.distinct
  }

  /**
   * Gets the files from the field.  The field must be a File, a FileExtension, or a Seq or Set of either.
   * @param field Field to get files.
   * @return for the field.
   */
  def getFieldFiles(field: ArgumentSource): Seq[File] = {
    var files: Seq[File] = Nil
    CollectionUtils.foreach(getFieldValue(field), (fieldValue) => {
      val file = fieldValueToFile(field, fieldValue)
      if (file != null)
        files :+= file
    })
    files.distinct
  }

  /**
   *  Gets the file from the field.  The field must be a File or a FileExtension and not a Seq or Set.
   * @param field Field to get the file.
   * @return for the field.
   */
  def getFieldFile(field: ArgumentSource): File =
    fieldValueToFile(field, getFieldValue(field))

  /**
   * Converts the field value to a file.  The field must be a File or a FileExtension.
   * @param field Field to get the file.
   * @param value Value of the File or FileExtension or null.
   * @return Null if value is null, otherwise the File.
   * @throws QException if the value is not a File or FileExtension.
   */
  private def fieldValueToFile(field: ArgumentSource, value: Any): File = value match {
    case file: File => file
    case null => null
    case unknown => throw new QException("Non-file found.  Try removing the annotation, change the annotation to @Argument, or extend File with FileExtension: %s: %s".format(field.field, unknown))
  }

  /**
   * After a function is frozen no more updates are allowed by the user.
   * The function is allow to make necessary updates internally to make sure
   * the inputs and outputs will be equal to other inputs and outputs.
   */
  final def freeze() {
    freezeFieldValues()
    canonFieldValues()
  }

  /**
   * Sets all field values.
   */
  def freezeFieldValues() {
    if (jobName == null)
      jobName = qSettings.runName + "-" + this.addOrder.mkString("-")

    if (jobOutputFile == null) {
      /*If the outputFile has been set to an absolute path, respect that.
        Otherwise, place it in (possibly a subdirectory of) the log directory
        The relative case is first as it's arguably the most common condition
      */
      jobOutputFile = firstOutput match {
        case file: File if !IOUtils.isSpecialFile(file) && !file.isAbsolute =>
            val logDir : File = if (file.getParentFile == null) qSettings.logDirectory else new File(qSettings.logDirectory, file.getParent)
            new File(logDir, file.getName + ".out")

        case file: File if !IOUtils.isSpecialFile(file) && file.isAbsolute =>
            new File(file.getParentFile, file.getName + ".out")

        case _ =>
          new File(qSettings.logDirectory, jobName + ".out")
      }
    }

    if (jobTempDir == null)
      jobTempDir = qSettings.tempDirectory

    if (jobLocalDir == null)
      jobLocalDir = jobTempDir

    if (jobPriority.isEmpty)
      jobPriority = qSettings.jobPriority

    // Do not set the temp and local dir relative to the command directory
    jobTempDir = IOUtils.absolute(jobTempDir)
    jobLocalDir = IOUtils.absolute(jobLocalDir)

    absoluteCommandDirectory()
  }

  /**
   * If the command directory is relative, insert the run directory ahead of it.
   */
  def absoluteCommandDirectory() {
    commandDirectory = IOUtils.absolute(qSettings.runDirectory, commandDirectory)
  }

  /**
   * Makes all field values canonical so that the graph can match the
   * inputs of one function to the output of another using equals().
   */
  def canonFieldValues() {
    for (field <- this.functionFields) {
      var fieldValue = this.getFieldValue(field)
      fieldValue = CollectionUtils.updated(fieldValue, canon).asInstanceOf[AnyRef]
      this.setFieldValue(field, fieldValue)
    }

    this.jobOutputFile = canon(this.jobOutputFile).asInstanceOf[File]
    if (this.jobErrorFile != null)
      this.jobErrorFile = canon(this.jobErrorFile).asInstanceOf[File]
  }

  /**
   * Set value to a uniform value across functions.
   * Base implementation changes any relative path to an absolute path.
   * @param value to be updated
   * @return the modified value, or a copy if the value is immutable
   */
  protected def canon(value: Any) = {
    value match {
      case file: File => IOUtils.absolute(commandDirectory, file)
      case x => x
    }
  }

  /**
   * Scala sugar type for checking annotation required and exclusiveOf.
   */
  private type ArgumentAnnotation = {
    def required(): Boolean
    def exclusiveOf(): String
    def doc(): String
  }

  /**
   * Returns the isRequired value from the field.
   * @param field Field to check.
   * @param annotation Annotation.
   * @return the isRequired value from the field annotation.
   */
  private def isRequired(field: ArgumentSource, annotation: Class[_ <: Annotation]) =
    ReflectionUtils.getAnnotation(field.field, annotation).asInstanceOf[ArgumentAnnotation].required()

  /**
   * Returns an array of ArgumentSources from functionFields listed in the exclusiveOf of the original field
   * @param field Field to check.
   * @param annotation Annotation.
   * @return the Array[ArgumentSource] that may be set instead of the field.
   */
  private def exclusiveOf(field: ArgumentSource, annotation: Class[_ <: Annotation]) =
    ReflectionUtils.getAnnotation(field.field, annotation).asInstanceOf[ArgumentAnnotation].exclusiveOf()
            .split(",").map(_.trim).filter(_.length > 0)
            .map(fieldName => functionFields.find(fieldName == _.field.getName) match {
      case Some(x) => x
      case None => throw new QException("Unable to find exclusion field %s on %s".format(fieldName, this.getClass.getSimpleName))
    })

  /**
   * Returns the doc value from the field.
   * @param field Field to check.
   * @param annotation Annotation.
   * @return the doc value from the field annotation.
   */
  private def doc(field: ArgumentSource, annotation: Class[_ <: Annotation]) =
    ReflectionUtils.getAnnotation(field.field, annotation).asInstanceOf[ArgumentAnnotation].doc()

  /**
   * Returns true if the field has a value.
   * @param source Field to check for a value.
   * @return true if the field has a value.
   */
  protected def hasFieldValue(source: ArgumentSource) = this.hasValue(this.getFieldValue(source))

  /**
   * Returns false if the value is null or an empty collection.
   * @param param Value to test for null, or a collection to test if it is empty.
   * @return false if the value is null, or false if the collection is empty, otherwise true.
   */
  protected def hasValue(param: Any) = CollectionUtils.isNotNullOrNotEmpty(param)

  /**
   * Gets the value of a field.
   * @param source Field to get the value for.
   * @return value of the field.
   */
  def getFieldValue(source: ArgumentSource) = ClassFieldCache.getFieldValue(this, source)

  /**
   * Gets the value of a field.
   * @param source Field to set the value for.
   * @return value of the field.
   */
  def setFieldValue(source: ArgumentSource, value: Any) = ClassFieldCache.setFieldValue(this, source, value)
}

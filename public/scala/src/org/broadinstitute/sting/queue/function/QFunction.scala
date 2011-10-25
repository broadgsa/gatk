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

package org.broadinstitute.sting.queue.function

import java.io.File
import java.lang.annotation.Annotation
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.{QException, QSettings}
import collection.JavaConversions._
import org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.utils.io.IOUtils

/**
 * The base interface for all functions in Queue.
 * Inputs and outputs are specified as Sets of values.
 * Inputs are matched to other outputs by using .equals()
 */
trait QFunction extends Logging with QJobReport {
  /** A short description of this step in the graph */
  var analysisName: String = "<function>"

  /** Prefix for automatic job name creation */
  var jobNamePrefix: String = _

  /** The name name of the job */
  var jobName: String = _

  /** Default settings */
  var qSettings: QSettings = _

  /** Directory to run the command in. */
  var commandDirectory: File = new File(".")

  /** Temporary directory to write any files */
  var jobTempDir: File = null

  /** Order the function was added to the graph. */
  var addOrder: List[Int] = Nil

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

  /**
   * If true and isIntermediate is true, the files listed
   * via outputs will deleted after the command completes.
   */
  var deleteIntermediateOutputs = true

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
    function.analysisName = this.analysisName
    function.jobName = this.jobName
    function.qSettings = this.qSettings
    function.commandDirectory = this.commandDirectory
    function.jobTempDir = this.jobTempDir
    function.addOrder = this.addOrder
    function.jobPriority = this.jobPriority
    function.jobRestartable = this.jobRestartable
    function.updateJobRun = this.updateJobRun
    function.isIntermediate = this.isIntermediate
    function.deleteIntermediateOutputs = this.deleteIntermediateOutputs
    function.reportGroup = this.reportGroup
    function.reportFeatures = this.reportFeatures
  }

  /** File to redirect any output.  Defaults to <jobName>.out */
  @Output(doc="File to redirect any output", required=false)
  @Gather(classOf[SimpleTextGatherFunction])
  var jobOutputFile: File = _

  /** File to redirect any errors.  Defaults to <jobName>.out */
  @Output(doc="File to redirect any errors", required=false)
  @Gather(classOf[SimpleTextGatherFunction])
  var jobErrorFile: File = _

  /**
   * Description of this command line function.
   */
  def description: String

  /**
   * The function description in .dot files
   */
  def dotString = jobName + " => " + description

  /**
   * Returns true if the function is done, false if it's
   * not done and None if the done status is unknown.
   */
  def isDone = {
    val files = doneOutputs
    if (files.size == 0)
      None
    else
      Some(files.forall(_.exists))
  }

  /**
   * Returns true if the function has failed, false if it
   * has not failed and None if the fail status is unknown.
   */
  def isFail = {
    val files = failOutputs
    if (files.size == 0)
      None
    else
      Some(files.exists(_.exists))
  }

  /**
   * Returns true if the file is a log file for this function.
   */
  protected def isLogFile(file: File) =
    file == jobOutputFile || file == jobErrorFile

  /**
   * Returns the output files for this function.
   * @return Set[File] outputs for this function.
   */
  private def statusPaths =
    commandOutputs.map(file => file.getParentFile + "/." + file.getName)
  
  /**
   * Returns the output files for this function.
   * @return Set[File] outputs for this function.
   */
  def doneOutputs = statusPaths.map(path => new File(path + ".done"))

  /**
   * Returns the output files for this function.
   * @return Set[File] outputs for this function.
   */
  def failOutputs = statusPaths.map(path => new File(path + ".fail"))

  /** The complete list of fields on this CommandLineFunction. */
  def functionFields = QFunction.classFields(this.functionFieldClass).functionFields
  /** The @Input fields on this CommandLineFunction. */
  def inputFields = QFunction.classFields(this.functionFieldClass).inputFields
  /** The @Output fields on this CommandLineFunction. */
  def outputFields = QFunction.classFields(this.functionFieldClass).outputFields
  /** The @Argument fields on this CommandLineFunction. */
  def argumentFields = QFunction.classFields(this.functionFieldClass).argumentFields

  /**
   * Returns the class that should be used for looking up fields.
   */
  protected def functionFieldClass = this.getClass

  /**
   * Returns the input files for this function.
   * @return Set[File] inputs for this function.
   */
  def inputs = getFieldFiles(inputFields)

  /**
   * Returns the output files for this function.
   * @return Set[File] outputs for this function.
   */
  def outputs = getFieldFiles(outputFields)

  /**
   * Returns the non-log outputs for this function.
   * @return the non-log outputs for this function.
   */
  def commandOutputs = outputs.filterNot(file => isLogFile(file))

  /**
   * Returns the set of directories where files may be written.
   */
  def outputDirectories = {
    var dirs = Set.empty[File]
    dirs += commandDirectory
    dirs += jobTempDir
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
    commandOutputs.foreach(file => IOUtils.tryDelete(file))
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
   * @return List[String] names of fields missing values.
   */
  def missingFields: List[String] = {
    val missingInputs = missingFields(inputFields, classOf[Input])
    val missingOutputs = missingFields(outputFields, classOf[Output])
    val missingArguments = missingFields(argumentFields, classOf[Argument])
    (missingInputs | missingOutputs | missingArguments).toList.sorted
  }

  /**
   * Returns fields that do not have values which are required.
   * @param sources Fields to check.
   * @param annotation Annotation.
   * @return Set[String] names of fields missing values.
   */
  private def missingFields(sources: List[ArgumentSource], annotation: Class[_ <: Annotation]): Set[String] = {
    var missing = Set.empty[String]
    for (source <- sources) {
      if (isRequired(source, annotation))
        if (!hasFieldValue(source))
          if (!exclusiveOf(source, annotation).exists(otherSource => hasFieldValue(otherSource)))
            missing += "@%s: %s - %s".format(annotation.getSimpleName, source.field.getName, doc(source, annotation))
    }
    missing
  }

  /**
   *  Gets the files from the fields.  The fields must be a File, a FileExtension, or a List or Set of either.
   * @param fields Fields to get files.
   * @return Set[File] for the fields.
   */
  private def getFieldFiles(fields: List[ArgumentSource]): Set[File] = {
    var files = Set.empty[File]
    for (field <- fields)
      files ++= getFieldFiles(field)
    files
  }

  /**
   * Gets the files from the field.  The field must be a File, a FileExtension, or a List or Set of either.
   * @param fields Field to get files.
   * @return Set[File] for the field.
   */
  def getFieldFiles(field: ArgumentSource): Set[File] = {
    var files = Set.empty[File]
    CollectionUtils.foreach(getFieldValue(field), (fieldValue) => {
      val file = fieldValueToFile(field, fieldValue)
      if (file != null)
        files += file
    })
    files
  }

  /**
   *  Gets the file from the field.  The field must be a File or a FileExtension and not a List or Set.
   * @param field Field to get the file.
   * @return File for the field.
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
    if (jobNamePrefix == null)
      jobNamePrefix = qSettings.jobNamePrefix

    if (jobName == null)
      jobName = QFunction.nextJobName(jobNamePrefix)

    if (jobOutputFile == null)
      jobOutputFile = new File(jobName + ".out")

    if (jobTempDir == null)
      jobTempDir = qSettings.tempDirectory

    if (jobPriority.isEmpty)
      jobPriority = qSettings.jobPriority

    // Do not set the temp dir relative to the command directory
    jobTempDir = IOUtils.absolute(jobTempDir)

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
   * @param value Value to test for null, or a collection to test if it is empty.
   * @return false if the value is null, or false if the collection is empty, otherwise true.
   */
  protected def hasValue(param: Any) = CollectionUtils.isNotNullOrNotEmpty(param)

  /**
   * Gets the value of a field.
   * @param source Field to get the value for.
   * @return value of the field.
   */
  def getFieldValue(source: ArgumentSource) = ReflectionUtils.getValue(invokeObj(source), source.field)

  /**
   * Gets the value of a field.
   * @param source Field to set the value for.
   * @return value of the field.
   */
  def setFieldValue(source: ArgumentSource, value: Any) = ReflectionUtils.setValue(invokeObj(source), source.field, value)

  /**
   * Walks gets the fields in this object or any collections in that object
   * recursively to find the object holding the field to be retrieved or set.
   * @param source Field find the invoke object for.
   * @return Object to invoke the field on.
   */
  private def invokeObj(source: ArgumentSource) = source.parentFields.foldLeft[AnyRef](this)(ReflectionUtils.getValue(_, _))
}

object QFunction {
  /** Job index counter for this run of Queue. */
  private var jobIndex = 0

  var parsingEngine: ParsingEngine = _

  /**
   * Returns the next job name using the prefix.
   * @param prefix Prefix of the job name.
   * @return the next job name.
   */
  private def nextJobName(prefix: String) = {
    jobIndex += 1
    prefix + "-" + jobIndex
  }

  /**
   * The list of fields defined on a class
   * @param clazz The class to lookup fields.
   */
  private class ClassFields(clazz: Class[_]) {
    /** The complete list of fields on this CommandLineFunction. */
    val functionFields: List[ArgumentSource] = parsingEngine.extractArgumentSources(clazz).toList
    /** The @Input fields on this CommandLineFunction. */
    val inputFields = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Input]))
    /** The @Output fields on this CommandLineFunction. */
    val outputFields = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Output]))
    /** The @Argument fields on this CommandLineFunction. */
    val argumentFields = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Argument]))
  }

  /**
   * The mapping from class to fields.
   */
  private var classFieldsMap = Map.empty[Class[_], ClassFields]

  /**
   * Returns the field on clazz.
   * @param clazz Class to search.
   * @param name Name of the field to return.
   * @return Argument source for the field.
   */
  def findField(clazz: Class[_], name: String) = {
    classFields(clazz).functionFields.find(_.field.getName == name) match {
      case Some(source) => source
      case None => throw new QException("Could not find a field on class %s with name %s".format(clazz, name))
    }
  }

  /**
   * Returns the fields for a class.
   * @param clazz Class to retrieve fields for.
   * @return the fields for the class.
   */
  private def classFields(clazz: Class[_]) = {
    classFieldsMap.get(clazz) match {
      case Some(classFields) => classFields
      case None =>
        val classFields = new ClassFields(clazz)
        classFieldsMap += clazz -> classFields
        classFields
    }
  }
}

package org.broadinstitute.sting.queue.function

import java.io.File
import java.lang.annotation.Annotation
import org.broadinstitute.sting.commandline._
import org.broadinstitute.sting.queue.util.{CollectionUtils, IOUtils, ReflectionUtils}
import org.broadinstitute.sting.queue.{QException, QSettings}
import collection.JavaConversions._

/**
 * The base interface for all functions in Queue.
 * Inputs and outputs are specified as Sets of values.
 * Inputs are matched to other outputs by using .equals()
 */
trait QFunction {
  /**
   * Analysis function name
   */
  var analysisName: String = _

  /** Default settings */
  var qSettings: QSettings = _

  /** Directory to run the command in. */
  var commandDirectory: File = IOUtils.CURRENT_DIR

  /**
   * Description of this command line function.
   */
  def description: String

  /**
   * The function description in .dot files
   */
  def dotString = ""

  protected def useStatusOutput(file: File): Boolean

  /**
   * Returns the output files for this function.
   * @return Set[File] outputs for this function.
   */
  private def statusPaths = outputs
          .filter(file => useStatusOutput(file))
          .map(file => file.getParentFile + "/." + file.getName)
  
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
  lazy val functionFields: List[ArgumentSource] = ParsingEngine.extractArgumentSources(this.getClass).toList
  /** The @Input fields on this CommandLineFunction. */
  lazy val inputFields = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Input]))
  /** The @Output fields on this CommandLineFunction. */
  lazy val outputFields = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Output]))
  /** The @Argument fields on this CommandLineFunction. */
  lazy val argumentFields = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Argument]))

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
   * Resets the field to the temporary directory.
   * @param field Field to get and set the file.
   * @param tempDir new root for the file.
   */
  def resetFieldFile(field: ArgumentSource, tempDir: File): File = {
    getFieldValue(field) match {
      case fileExtension: FileExtension => {
        val newFile = IOUtils.resetParent(tempDir, fileExtension)
        val newFileExtension = fileExtension.withPath(newFile.getPath)
        setFieldValue(field, newFileExtension)
        newFileExtension
      }
      case file: File => {
        if (file.getClass != classOf[File])
          throw new QException("Extensions of file must also extend with FileExtension so that the path can be modified.");
        val newFile = IOUtils.resetParent(tempDir, file)
        setFieldValue(field, newFile)
        newFile
      }
      case null => null
      case unknown =>
        throw new QException("Unable to set file from %s: %s".format(field, unknown))
    }
  }


  /**
   * After a function is frozen no more updates are allowed by the user.
   * The function is allow to make necessary updates internally to make sure
   * the inputs and outputs will be equal to other inputs and outputs.
   */
  final def freeze = {
    freezeFieldValues
    canonFieldValues
  }

  def freezeFieldValues = {
    commandDirectory = IOUtils.subDir(IOUtils.CURRENT_DIR, commandDirectory)
  }

  /**
   * Makes all field values canonical so that the graph can match the
   * inputs of one function to the output of another using equals().
   */
  def canonFieldValues = {
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
      case fileExtension: FileExtension =>
        val newFile = absolute(fileExtension);
        val newFileExtension = fileExtension.withPath(newFile.getPath)
        newFileExtension
      case file: File =>
        if (file.getClass != classOf[File])
          throw new QException("Extensions of file must also extend with FileExtension so that the path can be modified.");
        absolute(file)
      case x => x
    }
  }

  /**
   * Returns the absolute path to the file relative to the job command directory.
   * @param file File to root relative to the command directory if it is not already absolute.
   * @return The absolute path to file.
   */
  private def absolute(file: File) = IOUtils.subDir(commandDirectory, file)


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
    ReflectionUtils.getAnnotation(field.field, annotation).asInstanceOf[ArgumentAnnotation].required

  /**
   * Returns an array of ArgumentSources from functionFields listed in the exclusiveOf of the original field
   * @param field Field to check.
   * @param annotation Annotation.
   * @return the Array[ArgumentSource] that may be set instead of the field.
   */
  private def exclusiveOf(field: ArgumentSource, annotation: Class[_ <: Annotation]) =
    ReflectionUtils.getAnnotation(field.field, annotation).asInstanceOf[ArgumentAnnotation].exclusiveOf
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
    ReflectionUtils.getAnnotation(field.field, annotation).asInstanceOf[ArgumentAnnotation].doc

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

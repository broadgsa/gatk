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

package org.broadinstitute.gatk.queue.util

import org.broadinstitute.gatk.utils.commandline._
import scala.Some
import org.broadinstitute.gatk.queue.QException
import collection.JavaConversions._
import java.io.File

/**
 * Utilities and a static cache of argument fields for various classes populated by the parsingEngine.
 * Because this class works with the ParsingEngine it can walk @ArgumentCollection hierarchies.
 */
object ClassFieldCache {
  var parsingEngine: ParsingEngine = _


  //
  // Field caching
  //

  /**
   * The list of fields defined on a class
   * @param clazz The class to lookup fields.
   */
  private class ClassFields(clazz: Class[_]) {
    /** The complete list of fields on this CommandLineFunction. */
    val functionFields: Seq[ArgumentSource] = parsingEngine.extractArgumentSources(clazz).toSeq
    /** The @Input fields on this CommandLineFunction. */
    val inputFields: Seq[ArgumentSource] = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Input]))
    /** The @Output fields on this CommandLineFunction. */
    val outputFields: Seq[ArgumentSource] = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Output]))
    /** The @Argument fields on this CommandLineFunction. */
    val argumentFields: Seq[ArgumentSource] = functionFields.filter(source => ReflectionUtils.hasAnnotation(source.field, classOf[Argument]))
  }

  /**
   * The mapping from class to fields.
   */
  private var classFieldsMap = Map.empty[Class[_], ClassFields]

  /**
   * Returns the fields for a class.
   * @param clazz Class to retrieve fields for.
   * @return the fields for the class.
   */
  private def classFields(clazz: Class[_]): ClassFields = {
    classFieldsMap.get(clazz) match {
      case Some(classFields) => classFields
      case None =>
        val classFields = new ClassFields(clazz)
        classFieldsMap += clazz -> classFields
        classFields
    }
  }

  /**
   * Returns the field on clazz.
   * @param clazz Class to search.
   * @param name Name of the field to return.
   * @return Argument source for the field.
   */
  def findField(clazz: Class[_], name: String): ArgumentSource = {
    classFields(clazz).functionFields.find(_.field.getName == name) match {
      case Some(source) => source
      case None => throw new QException("Could not find a field on class %s with name %s".format(clazz, name))
    }
  }

  /**
   * Returns the Seq of fields for a QFunction class.
   * @param clazz Class to retrieve fields for.
   * @return the fields of the class.
   */
  def classFunctionFields(clazz: Class[_]): Seq[ArgumentSource] = classFields(clazz).functionFields

  /**
   * Returns the Seq of inputs for a QFunction class.
   * @param clazz Class to retrieve inputs for.
   * @return the inputs of the class.
   */
  def classInputFields(clazz: Class[_]): Seq[ArgumentSource] = classFields(clazz).inputFields

  /**
   * Returns the Seq of outputs for a QFunction class.
   * @param clazz Class to retrieve outputs for.
   * @return the outputs of the class.
   */
  def classOutputFields(clazz: Class[_]): Seq[ArgumentSource] = classFields(clazz).outputFields

  /**
   * Returns the Seq of arguments for a QFunction class.
   * @param clazz Class to retrieve arguments for.
   * @return the arguments of the class.
   */
  def classArgumentFields(clazz: Class[_]): Seq[ArgumentSource] = classFields(clazz).argumentFields


  //
  // get/set fields as AnyRef
  //

  /**
   * Gets the value of a field.
   * @param obj Top level object storing the source info.
   * @param source Field to get the value for.
   * @return value of the field.
   */
  def getFieldValue(obj: AnyRef, source: ArgumentSource) = ReflectionUtils.getValue(invokeObj(obj, source), source.field)

  /**
   * Gets the value of a field.
   * @param obj Top level object storing the source info.
   * @param source Field to set the value for.
   * @return value of the field.
   */
  def setFieldValue(obj: AnyRef, source: ArgumentSource, value: Any) = ReflectionUtils.setValue(invokeObj(obj, source), source.field, value)

  /**
   * Walks gets the fields in this object or any collections in that object
   * recursively to find the object holding the field to be retrieved or set.
   * @param obj Top level object storing the source info.
   * @param source Field find the invoke object for.
   * @return Object to invoke the field on.
   */
  private def invokeObj(obj: AnyRef, source: ArgumentSource) = source.parentFields.foldLeft[AnyRef](obj)(ReflectionUtils.getValue(_, _))


  //
  // get/set fields as java.io.File
  //

  /**
   *  Gets the files from the fields.  The fields must be a File, a FileExtension, or a Seq or Set of either.
   * @param obj Top level object storing the source info.
   * @param fields Fields to get files.
   * @return for the fields.
   */
  def getFieldFiles(obj: AnyRef, fields: Seq[ArgumentSource]): Seq[File] = {
    var files: Seq[File] = Nil
    for (field <- fields)
      files ++= getFieldFiles(obj, field)
    files.distinct
  }

  /**
   * Gets the files from the field.  The field must be a File, a FileExtension, or a Seq or Set of either.
   * @param obj Top level object storing the source info.
   * @param field Field to get files.
   * @return for the field.
   */
  def getFieldFiles(obj: AnyRef, field: ArgumentSource): Seq[File] = {
    var files: Seq[File] = Nil
    CollectionUtils.foreach(getFieldValue(obj, field), (fieldValue) => {
      val file = fieldValueToFile(field, fieldValue)
      if (file != null)
        files :+= file
    })
    files.distinct
  }

  /**
   *  Gets the file from the field.  The field must be a File or a FileExtension and not a Seq or Set.
   * @param obj Top level object storing the source info.
   * @param field Field to get the file.
   * @return for the field.
   */
  def getFieldFile(obj: AnyRef, field: ArgumentSource): File =
    fieldValueToFile(field, getFieldValue(obj, field))

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


  //
  // other utilities
  //

  /**
   * Retrieves the fullName of the argument
   * @param field ArgumentSource to check
   * @return Full name of the argument source
   */
  def fullName(field: ArgumentSource) = field.createArgumentDefinitions().get(0).fullName
}

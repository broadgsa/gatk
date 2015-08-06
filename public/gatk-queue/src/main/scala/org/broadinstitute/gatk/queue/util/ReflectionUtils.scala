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

import org.broadinstitute.gatk.queue.QException
import java.lang.annotation.Annotation
import java.lang.reflect.{ParameterizedType, Field}
import org.broadinstitute.gatk.utils.commandline.ClassType
import org.broadinstitute.gatk.utils.classloader.JVMUtils

/**
 * A collection of scala extensions to the GATK JVMUtils.
 */
object ReflectionUtils {

  /**
   * Returns true if field has the annotation.
   * @param field Field to check.
   * @param annotation Class of the annotation to look for.
   * @return true if field has the annotation.
   */
  def hasAnnotation(field: Field, annotation: Class[_ <: Annotation]) = field.getAnnotation(annotation) != null

  /**
   * Returns true if clazz or one of its superclasses has the annotation.
   * @param clazz Class to check.
   * @param annotation Class of the annotation to look for.
   * @return true if field has the annotation.
   */
  def hasAnnotation(clazz: Class[_], annotation: Class[_ <: Annotation]) = {
    var foundAnnotation = false
    while (!foundAnnotation && clazz != null)
      foundAnnotation = (clazz.getAnnotation(annotation) != null)
    foundAnnotation
  }

  /**
   * Gets the annotation or throws an exception if the annotation is not found.
   * @param field Field to check.
   * @param annotation Class of the annotation to look for.
   * @return The annotation.
   */
  def getAnnotation[T <: Annotation](field: Field, annotation: Class[T]): T = {
    field.getAnnotation(annotation) match {
      case null =>
        throw new QException("Field %s is missing annotation %s".format(field, annotation))
      case fieldAnnotation => fieldAnnotation.asInstanceOf[T]
    }
  }

  /**
   * Gets the annotation or throws an exception if the annotation is not found.
   * @param clazz Class to check.
   * @param annotation Class of the annotation to look for.
   * @return The annotation.
   */
  def getAnnotation[T <: Annotation](clazz: Class[_], annotation: Class[T]): T = {
    var result: T = null.asInstanceOf[T]
    while (result == null && clazz != null)
      result = clazz.getAnnotation(annotation)
    if (result == null)
      throw new QException("Class %s is missing annotation %s".format(clazz, annotation))
    result
  }

  /**
   * Returns all the declared fields on a class in order of sub type to super type.
   * @param clazz Base class to start looking for fields.
   * @return Seq[Field] found on the class and all super classes.
   */
  def getAllFields(clazz: Class[_]) = getAllTypes(clazz).map(_.getDeclaredFields).flatMap(_.toSeq)

  /**
   * Gets all the types on a class in order of sub type to super type.
   * @param clazz Base class.
   * @return Seq[Class] including the class and all super classes.
   */
  def getAllTypes(clazz: Class[_]) = {
    var types = Seq.empty[Class[_]]
    var c = clazz
    while (c != null) {
      types :+= c
      c = c.getSuperclass
    }
    types
  }

  /**
   * Gets a field value using reflection.
   * Attempts to use the scala getter then falls back to directly accessing the field.
   * @param obj Object to inspect.
   * @param field Field to retrieve.
   * @return The field value.
   */
  def getValue(obj: AnyRef, field: Field): AnyRef =
    try {
      field.getDeclaringClass.getMethod(field.getName).invoke(obj)
    } catch {
      case e: NoSuchMethodException => JVMUtils.getFieldValue(field, obj)
    }

  /**
   * Sets a field value using reflection.
   * Attempts to use the scala setter then falls back to directly accessing the field.
   * @param obj Object to inspect.
   * @param field Field to set.
   * @param value The new field value.
   */
  def setValue(obj: AnyRef, field: Field, value: Any) =
    try {
      field.getDeclaringClass.getMethod(field.getName+"_$eq", field.getType).invoke(obj, value.asInstanceOf[AnyRef])
    } catch {
      case e: NoSuchMethodException => JVMUtils.setFieldValue(field, obj, value)
    }

  /**
   * Returns the collection type of a field or throws an exception if the field contains more than one parameterized type, or the collection type cannot be found.
   * @param field Field to retrieve the collection type.
   * @return The collection type for the field.
   */
  def getCollectionType(field: Field) = {
    getGenericTypes(field) match {
      case Some(classes) =>
        if (classes.length > 1)
          throw new IllegalArgumentException("Field contains more than one generic type: " + field)
        classes(0)
      case None =>
        throw new QException("Generic type not set for collection.  Did it declare an @ClassType?: " + field)
    }
  }

  /**
   * Returns the generic types for a field or None.
   * @param field Field to retrieve the collection type.
   * @return The array of classes that are in the collection type, or None if the type cannot be found.
   */
  private def getGenericTypes(field: Field): Option[Array[Class[_]]] = {
    // TODO: Refactor: based on java code in org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor
    // If this is a parameterized collection, find the contained type.  If blow up if only one type exists.
    if (hasAnnotation(field, classOf[ClassType])) {
      Some(Array(getAnnotation(field, classOf[ClassType]).value))
    } else if (field.getGenericType.isInstanceOf[ParameterizedType]) {
      val parameterizedType = field.getGenericType.asInstanceOf[ParameterizedType]
      Some(parameterizedType.getActualTypeArguments.map(_.asInstanceOf[Class[_]]))
    } else None
  }
}

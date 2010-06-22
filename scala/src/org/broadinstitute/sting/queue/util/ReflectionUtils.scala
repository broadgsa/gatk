package org.broadinstitute.sting.queue.util

import org.broadinstitute.sting.queue.QException
import java.lang.annotation.Annotation
import scala.concurrent.JavaConversions
import scala.concurrent.JavaConversions._
import scala.collection.immutable.ListMap
import java.lang.reflect.{ParameterizedType, Field}

object ReflectionUtils {
  def hasAnnotation(field: Field, annotation: Class[_ <: Annotation]) = field.getAnnotation(annotation) != null

  def getAnnotation[T <: Annotation](field: Field, annotation: Class[T]): T = {
    if (!hasAnnotation(field, annotation))
      throw new QException("Field %s is missing annotation %s".format(field, annotation))
    field.getAnnotation(annotation).asInstanceOf[T]
  }
  
  def getAllFields(clazz: Class[_]) = getAllTypes(clazz).map(_.getDeclaredFields).flatMap(_.toList)

  def filterFields(fields: List[Field], annotation: Class[_ <: Annotation]) = fields.filter(field => hasAnnotation(field, annotation))

  def getFieldNamesValues(obj: AnyRef, fields: List[Field]) =
    ListMap(fields.map(field => (field.getName -> fieldGetter(field).invoke(obj))) :_*)

  def getAllTypes(clazz: Class[_]) = {
    var types = List.empty[Class[_]]
      var c = clazz
      while (c != null) {
        types :+= c
        c = c.getSuperclass
      }
    types
  }

  def getValue(obj: AnyRef, field: Field) = fieldGetter(field).invoke(obj)
  def setValue(obj: AnyRef, field: Field, value: Any) = fieldSetter(field).invoke(obj, value.asInstanceOf[AnyRef])

  def addOrUpdateWithStringValue(obj: AnyRef, field: Field, value: String) = {
    val getter = fieldGetter(field)
    val setter = fieldSetter(field)

    if (classOf[Seq[_]].isAssignableFrom(field.getType)) {

      val fieldType = getCollectionType(field)
      val typeValue = coerce(fieldType, value)

      var list = getter.invoke(obj).asInstanceOf[Seq[_]]
      list :+= typeValue
      setter.invoke(obj, list)

    } else if (classOf[Option[_]].isAssignableFrom(field.getType)) {

      val fieldType = getCollectionType(field)
      val typeValue = coerce(fieldType, value)

      setter.invoke(obj, Some(typeValue))

    } else {

      val fieldType = field.getType
      val typeValue = coerce(fieldType, value)

      setter.invoke(obj, typeValue.asInstanceOf[AnyRef])
    }
  }

  private def getCollectionType(field: Field) = {
    getGenericTypes(field) match {
      case Some(classes) =>
        if (classes.length > 1)
          throw new IllegalArgumentException("Field contains more than one generic type: " + field)
        classes(0)
      case None =>
        if (!field.isAnnotationPresent(classOf[ClassType]))
          throw new QException("@ClassType must be specified for unparameterized field: " + field)
        field.getAnnotation(classOf[ClassType]).asInstanceOf[ClassType].value
    }
  }

  private def getGenericTypes(field: Field) = {
    // TODO: Refactor: based on java code in org.broadinstitute.sting.commandline.ArgumentTypeDescriptor
    // If this is a parameterized collection, find the contained type.  If blow up if only one type exists.
    if (field.getGenericType.isInstanceOf[ParameterizedType]) {
      val parameterizedType = field.getGenericType.asInstanceOf[ParameterizedType]
      Some(parameterizedType.getActualTypeArguments.map(_.asInstanceOf[Class[_]]))
    }
    else None
  }

  private def fieldGetter(field: Field) =
    try {
      field.getDeclaringClass.getMethod(field.getName)
    } catch {
      case e: NoSuchMethodException => throw new QException("Field may be private?  Unable to find getter for field: " + field)
    }

  private def fieldSetter(field: Field) =
    try {
      field.getDeclaringClass.getMethod(field.getName+"_$eq", field.getType)
    } catch {
      case e: NoSuchMethodException => throw new QException("Field may be a val instead of var?  Unable to find setter for field: " + field)
    }

  private def coerce(clazz: Class[_], value: String) = {
    if (classOf[String] == clazz) value
    else if (classOf[Boolean] == clazz) value.toBoolean
    else if (classOf[Byte] == clazz) value.toByte
    else if (classOf[Short] == clazz) value.toShort
    else if (classOf[Int] == clazz) value.toInt
    else if (classOf[Long] == clazz) value.toLong
    else if (classOf[Float] == clazz) value.toFloat
    else if (classOf[Double] == clazz) value.toDouble
    else if (hasStringConstructor(clazz))
      clazz.getConstructor(classOf[String]).newInstance(value)
    else throw new QException("Unable to coerce value '%s' to type '%s'.".format(value, clazz))
  }

  private def hasStringConstructor(clazz: Class[_]) = {
    clazz.getConstructors.exists(constructor => {
      val parameters = constructor.getParameterTypes
      parameters.size == 1 && parameters.head == classOf[String]
    })
  }
}

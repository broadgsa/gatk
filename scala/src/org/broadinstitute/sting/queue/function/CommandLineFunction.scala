package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.util._
import java.lang.reflect.Field
import java.lang.annotation.Annotation
import org.broadinstitute.sting.commandline.{Input, Output, ArgumentDescription}

trait CommandLineFunction extends InputOutputFunction with DispatchFunction {
  var properties = Map.empty[String, String]

  def inputFieldsWithValues = inputFields.filter(hasFieldValue(_))
  def outputFieldsWithValues = outputFields.filter(hasFieldValue(_))

  /**
   * Sets parameters from the arg map.
   */
  override def freeze = {
    for ((name, value) <- properties) addOrUpdateWithStringValue(name, value)
    super.freeze
  }

  /**
   * Repeats parameters with a prefix/suffix if they are set otherwise returns "".
   * Skips null, Nil, None.  Unwraps Some(x) to x.  Everything else is called with x.toString.
   */
  protected def repeat(prefix: String, params: Seq[_], suffix: String = "", separator: String = "") =
    params.filter(param => hasValue(param)).map(param => prefix + toValue(param) + suffix).mkString(separator)

  /**
   * Returns parameter with a prefix/suffix if it is set otherwise returns "".
   * Does not output null, Nil, None.  Unwraps Some(x) to x.  Everything else is called with x.toString.
   */
  protected def optional(prefix: String, param: Any, suffix: String = "") =
    if (hasValue(param)) prefix + toValue(param) + suffix else ""

  def missingValues = {
    val missingInputs = missingFields(inputFields, classOf[Input])
    val missingOutputs = missingFields(outputFields, classOf[Output])
    missingInputs | missingOutputs
  }

  private def missingFields(fields: List[Field], annotation: Class[_ <: Annotation]) = {
    var missing = Set.empty[String]
    for (field <- fields) {
      if (isRequired(field, annotation))
        if (!hasValue(ReflectionUtils.getValue(this, field)))
          missing += field.getName
    }
    missing
  }

  private def isRequired(field: Field, annotation: Class[_ <: Annotation]) =
    new ArgumentDescription(field.getAnnotation(annotation)).required

  protected def hasFieldValue(field: Field) = hasValue(this.getFieldValue(field))

  private def hasValue(param: Any) = param match {
    case null => false
    case Nil => false
    case None => false
    case _ => true
  }

  private def toValue(param: Any): String = param match {
    case null => ""
    case Nil => ""
    case None => ""
    case Some(x) => x.toString
    case x => x.toString
  }
}

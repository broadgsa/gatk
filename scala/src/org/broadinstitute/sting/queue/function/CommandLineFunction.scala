package org.broadinstitute.sting.queue.function

import java.io.File
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.engine.{CommandLineRunner, QGraph}
import java.lang.reflect.Field

trait CommandLineFunction extends InputOutputFunction with DispatchFunction {
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
    val missingInputs = missingFields(inputFields)
    val missingOutputs = missingFields(outputFields)
    missingInputs | missingOutputs
  }

  private def missingFields(fields: List[Field]) = {
    var missing = Set.empty[String]
    for (field <- fields) {
      val isOptional = ReflectionUtils.hasAnnotation(field, classOf[Optional])
      if (!isOptional)
        if (!hasValue(ReflectionUtils.getValue(this, field)))
          missing += field.getName
    }
    missing
  }

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

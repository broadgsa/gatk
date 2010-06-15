package org.broadinstitute.sting.queue.function

import java.io.File
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.queue.engine.{CommandLineRunner, QGraph}

trait CommandLineFunction extends QFunction with DispatchFunction {

  /**
   * The command line to run locally or via grid computing.
   */
  def commandLine: String

  /**
   * The directory where the command should run.
   */
  @Internal
  var commandDirectory: File = new File(".")

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

  /**
   * Sets a field value using the name of the field.
   * Field must be annotated with @Input, @Output, or @Internal
   * @returns true if the value was found and set
   */
  def setValue(name: String, value: String) = {
    ReflectionUtils.getField(this, name) match {
      case Some(field) =>
        val isInput = ReflectionUtils.hasAnnotation(field, classOf[Input])
        val isOutput = ReflectionUtils.hasAnnotation(field, classOf[Output])
        val isInternal = ReflectionUtils.hasAnnotation(field, classOf[Internal])
        if (isInput || isOutput || isInternal) {
          ReflectionUtils.setValue(this, field, value)
        }
        true
      case None => false
    }
  }

  private lazy val fields = ReflectionUtils.getAllFields(this.getClass)
  private def internals = ReflectionUtils.getFieldsAnnotatedWith(this, fields, classOf[Internal])
  def inputs = ReflectionUtils.getFieldsAnnotatedWith(this, fields, classOf[Input])
  def outputs = ReflectionUtils.getFieldsAnnotatedWith(this, fields, classOf[Output])

  override def missingValues = {
    val missingInputs = missingFields(inputs)
    val missingOutputs = missingFields(outputs)
    missingInputs | missingOutputs
  }

  private def missingFields(fields: Map[String, Any]) = {
    var missing = Set.empty[String]
    for ((name, value) <- fields) {
      val isOptional = ReflectionUtils.getField(this, name) match {
        case Some(field) => ReflectionUtils.hasAnnotation(field, classOf[Optional])
        case None => false
      }
      if (!isOptional)
        if (!hasValue(value))
          missing += name
    }
    missing
  }

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

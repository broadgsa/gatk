package org.broadinstitute.sting.queue.function

import java.lang.reflect.Field
import org.broadinstitute.sting.queue.util._
import org.broadinstitute.sting.commandline.{Input, Output}

/**
 * A function with @Inputs and @Outputs tagging fields that can be set by the user in a QScript
 */
trait InputOutputFunction extends QFunction with Cloneable {
  def getFieldValue(field: Field) = ReflectionUtils.getValue(this, field)
  def setFieldValue(field: Field, value: Any) = ReflectionUtils.setValue(this, field, value)

  def functionFields: List[Field] = inputFields ::: outputFields
  def inputFields = ReflectionUtils.filterFields(fields, classOf[Input])
  def outputFields = ReflectionUtils.filterFields(fields, classOf[Output])

  private lazy val fields = ReflectionUtils.getAllFields(this.getClass)
  // TODO: Need to handle argument collections where field is not on THIS
  def inputs = CollectionUtils.removeNullOrEmpty(ReflectionUtils.getFieldValues(this, inputFields)).toSet
  def outputs = CollectionUtils.removeNullOrEmpty(ReflectionUtils.getFieldValues(this, outputFields)).toSet

  /**
   * Sets a field value using the name of the field.
   * Field must be annotated with @Input, @Output, or @Internal
   * @return true if the value was found and set
   */
  def addOrUpdateWithStringValue(name: String, value: String) = {
    fields.find(_.getName == name) match {
      case Some(field) =>
        val isInput = ReflectionUtils.hasAnnotation(field, classOf[Input])
        val isOutput = ReflectionUtils.hasAnnotation(field, classOf[Output])
        if (isInput || isOutput) {
          ReflectionUtils.addOrUpdateWithStringValue(this, field, value)
        }
        true
      // TODO: Need to handle argument collections where field is not on THIS
      case None => false
    }
  }

  def cloneFunction() = clone.asInstanceOf[this.type]
  // explicitly overriden so that trait function cloneFunction can use this.clone
  override protected def clone = super.clone

  /**
   * As the function is frozen, changes all fields to their canonical forms.
   */
  override def freeze = {
    for (field <- this.functionFields)
      mapField(field, canon)
    super.freeze
  }

  def mapField(field: Field, f: Any => Any): Any = {
    var fieldValue = this.getFieldValue(field)
    fieldValue = CollectionUtils.updated(fieldValue, f).asInstanceOf[AnyRef]
    this.setFieldValue(field, fieldValue)
    fieldValue
  }

  /**
   * Set value to a uniform value across functions.
   * The biggest example is file paths relative to the command directory in DispatchFunction
   */
  protected def canon(value: Any): Any = value
}
